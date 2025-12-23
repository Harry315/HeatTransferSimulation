# -*- coding: utf-8 -*-
"""
Full pipeline: STP → mesh/inp → heat analysis → HFL → keff

Steps:
1. From STP: generate 1/8 unit-cell geometry → mirror to full cell → build array → mesh → export inp
2. From inp: import model, convert to heat-transfer model, apply BCs & material, solve
3. From odb: read HFL on a given node set, compute average heat flux q
4. From q, L, ΔT: compute effective thermal conductivity keff
5. Clean intermediate files (keep odb and result csv)

Notes:
- Intended to run in `abaqus python` nogui mode (HPC / cluster scenarios)
- If you want GUI-based HFL extraction, you can additionally add `session`-based functions
"""

from abaqus import *
from abaqusConstants import *
from odbAccess import openOdb

import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior

import os
import time
import numpy as np
from caeModules import *
from driverUtils import executeOnCaeStartup


# ============================================================
# 1. STP → 1/8 unit cell → mirror → array → mesh → export inp
# ============================================================

def import_stp_to_part(stp_path, model_name, part_name):
    """Read STP, create a Part, and return max coordinates and mid coordinates"""
    step_file = mdb.openStep(stp_path, scaleFromFile=OFF)
    mdb.models[model_name].PartFromGeometryFile(
        name=part_name,
        geometryFile=step_file,
        combine=True,
        mergeSolidRegions=True,
        dimensionality=THREE_D,
        type=DEFORMABLE_BODY
    )

    p = mdb.models[model_name].parts[part_name]
    x_coords = [v.pointOn[0][0] for v in p.vertices]
    y_coords = [v.pointOn[0][1] for v in p.vertices]
    z_coords = [v.pointOn[0][2] for v in p.vertices]

    x_min, x_max = min(x_coords), max(x_coords)
    y_min, y_max = min(y_coords), max(y_coords)
    z_min, z_max = min(z_coords), max(z_coords)

    x_mid = 0.5 * (x_max + x_min)
    y_mid = 0.5 * (y_max + y_min)
    z_mid = 0.5 * (z_max + z_min)

    return (x_max, y_max, z_max, x_mid, y_mid, z_mid)


def build_box_for_eighth(model_name, x_max, x_center, y_max, y_center, z_max, z_center):
    """Create a box covering 1/8 of the space, used to cut a 1/8 unit cell from the full structure"""
    m = mdb.models[model_name]
    s = m.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(x_center, y_center), point2=(x_max, y_max))

    p = m.Part(name='box', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s, depth=(z_max - z_center))

    s.unsetPrimaryObject()
    del m.sketches['__profile__']


def cut_to_eighth(model_name, part_name):
    """Use the box and original structure to perform boolean cuts and obtain Part-8th"""
    m = mdb.models[model_name]
    a = m.rootAssembly
    a.DatumCsysByDefault(CARTESIAN)

    p_box = m.parts['box']
    a.Instance(name='box-1', part=p_box, dependent=ON)

    p_orig = m.parts[part_name]
    inst_orig_name = part_name + '-1'
    a.Instance(name=inst_orig_name, part=p_orig, dependent=ON)

    a.InstanceFromBooleanCut(
        name='Part-boxSub8th',
        instanceToBeCut=a.instances['box-1'],
        cuttingInstances=(a.instances[inst_orig_name],),
        originalInstances=SUPPRESS
    )

    a = m.rootAssembly
    a.Instance(name='box-2', part=p_box, dependent=ON)
    a.InstanceFromBooleanCut(
        name='Part-8th',
        instanceToBeCut=a.instances['box-2'],
        cuttingInstances=(a.instances['Part-boxSub8th-1'],),
        originalInstances=SUPPRESS
    )


def mesh_eighth(model_name, part_name='Part-8th',
                seed_size=0.4, dev_factor=0.08, min_size_factor=0.01):
    """Mesh Part-8th"""
    p = mdb.models[model_name].parts[part_name]
    c = p.cells
    picked_cells = c.getSequenceFromMask(mask=('[#1 ]',),)

    p.setMeshControls(regions=picked_cells, elemShape=TET, technique=FREE)

    elemType1 = mesh.ElemType(elemCode=C3D20R)
    elemType2 = mesh.ElemType(elemCode=C3D15)
    elemType3 = mesh.ElemType(elemCode=C3D10)
    p.setElementType(
        regions=(picked_cells,),
        elemTypes=(elemType1, elemType2, elemType3)
    )

    p.seedPart(size=seed_size, deviationFactor=dev_factor, minSizeFactor=min_size_factor)
    p.generateMesh()


def mirror_full_cell(model_name, x_center, y_center, z_center):
    """1/8 → 1/4 → 1/2 → full unit cell: mirror-finished"""
    m = mdb.models[model_name]
    a = m.rootAssembly

    # 1/8 + XY mirror
    p_8th = m.parts['Part-8th']
    p_8th_mxy = m.Part(
        name='Part-8th-MirrorXY',
        objectToCopy=p_8th,
        compressFeatureList=ON,
        mirrorPlane=XYPLANE
    )
    a.Instance(name='Part-8th-1', part=p_8th, dependent=ON)
    a.Instance(name='Part-8th-MirrorXY-1', part=p_8th_mxy, dependent=ON)
    a.InstanceFromBooleanMerge(
        name='Part-8th-copy',
        instances=(a.instances['Part-8th-1'], a.instances['Part-8th-MirrorXY-1']),
        mergeNodes=BOUNDARY_ONLY,
        nodeMergingTolerance=1e-06,
        domain=MESH,
        originalInstances=SUPPRESS
    )

    # 1/4 + YZ mirror
    p_8th_copy = m.parts['Part-8th-copy']
    p_8th_copy_myz = m.Part(
        name='Part-8th-copy-MirrorYZ',
        objectToCopy=p_8th_copy,
        compressFeatureList=ON,
        mirrorPlane=YZPLANE
    )
    a.Instance(name='Part-8th-copy-1', part=p_8th_copy, dependent=ON)
    a.Instance(name='Part-8th-copy-MirrorYZ-1', part=p_8th_copy_myz, dependent=ON)
    a.translate(instanceList=('Part-8th-copy-1',),
                vector=(-x_center, -y_center, -z_center))
    a.translate(instanceList=('Part-8th-copy-MirrorYZ-1',),
                vector=(+x_center, -y_center, -z_center))
    a.InstanceFromBooleanMerge(
        name='Part-4th',
        instances=(a.instances['Part-8th-copy-1'], a.instances['Part-8th-copy-MirrorYZ-1']),
        mergeNodes=BOUNDARY_ONLY,
        nodeMergingTolerance=1e-06,
        domain=MESH,
        originalInstances=SUPPRESS
    )

    # 1/2 + XZ mirror
    p_4th = m.parts['Part-4th']
    p_4th_mxz = m.Part(
        name='Part-4th-copy-MirrorXZ',
        objectToCopy=p_4th,
        compressFeatureList=ON,
        mirrorPlane=XZPLANE
    )
    a.Instance(name='Part-4th-1', part=p_4th, dependent=ON)
    a.Instance(name='Part-4th-copy-MirrorXZ-1', part=p_4th_mxz, dependent=ON)
    a.InstanceFromBooleanMerge(
        name='Part-2nd',
        instances=(a.instances['Part-4th-1'], a.instances['Part-4th-copy-MirrorXZ-1']),
        mergeNodes=BOUNDARY_ONLY,
        nodeMergingTolerance=1e-06,
        domain=MESH,
        originalInstances=SUPPRESS
    )

    # Full cell + XY mirror
    p_2nd = m.parts['Part-2nd']
    p_2nd_mxy = m.Part(
        name='Part-2nd-copy-MirrorXY',
        objectToCopy=p_2nd,
        compressFeatureList=ON,
        mirrorPlane=XYPLANE
    )
    a.Instance(name='Part-2nd-1', part=p_2nd, dependent=ON)
    a.Instance(name='Part-2nd-copy-MirrorXY-1', part=p_2nd_mxy, dependent=ON)
    a.InstanceFromBooleanMerge(
        name='mirror-finished',
        instances=(a.instances['Part-2nd-1'], a.instances['Part-2nd-copy-MirrorXY-1']),
        mergeNodes=BOUNDARY_ONLY,
        nodeMergingTolerance=1e-06,
        domain=MESH,
        originalInstances=SUPPRESS
    )


def build_array(model_name, arrayX_num, arrayY_num, arrayZ_num,
                x_len, y_len, z_len,
                x_center, y_center, z_center):
    """Build a 3D array of mirror-finished and merge into finalInstance"""
    m = mdb.models[model_name]
    a = m.rootAssembly

    # XY array
    a.LinearInstancePattern(
        instanceList=('mirror-finished-1',),
        direction1=(1.0, 0.0, 0.0),
        direction2=(0.0, 1.0, 0.0),
        number1=arrayX_num,
        number2=arrayY_num,
        spacing1=x_len,
        spacing2=y_len
    )

    base_inst = 'mirror-finished'
    inst_tuple = tuple(
        inst for inst in a.instances.values()
        if inst.name.startswith(base_inst)
    )
    a.InstanceFromBooleanMerge(
        name='array-finishedxy',
        instances=inst_tuple,
        mergeNodes=BOUNDARY_ONLY,
        nodeMergingTolerance=1e-06,
        domain=MESH,
        originalInstances=SUPPRESS
    )

    # Z array
    a.LinearInstancePattern(
        instanceList=('array-finishedxy-1',),
        direction1=(1.0, 0.0, 0.0),
        direction2=(0.0, 0.0, 1.0),
        number1=1,
        number2=arrayZ_num,
        spacing1=x_len,
        spacing2=z_len
    )

    base_inst = 'array-finishedxy-1'
    inst_tuple = tuple(
        inst for inst in a.instances.values()
        if inst.name.startswith(base_inst)
    )
    a.InstanceFromBooleanMerge(
        name='finalInstance',
        instances=inst_tuple,
        mergeNodes=BOUNDARY_ONLY,
        nodeMergingTolerance=1e-06,
        domain=MESH,
        originalInstances=SUPPRESS
    )

    # Translate to be centered
    a.translate(
        instanceList=('finalInstance-1',),
        vector=(
            x_center - arrayX_num * x_len * 0.5,
            y_center - arrayY_num * y_len * 0.5,
            z_center - arrayZ_num * z_len * 0.5
        )
    )


def clean_intermediate_instances(model_name):
    """Delete intermediate instances and keep only finalInstance-1"""
    a = mdb.models[model_name].rootAssembly
    a.regenerate()

    to_delete = tuple(
        inst for inst in a.instances.values()
        if (inst.name.startswith('mirror')
            or inst.name.lower().startswith('part')
            or inst.name.lower().startswith('array'))
        and not inst.name.startswith('finalInstance')
    )
    for inst in to_delete:
        a.deleteFeatures((inst.name,))


def write_inp_from_model(model_name, job_name):
    """Create a Job and write an inp file (no solving)"""
    if job_name in mdb.jobs.keys():
        del mdb.jobs[job_name]

    mdb.Job(
        name=job_name,
        model=model_name,
        description='mesh export job',
        type=ANALYSIS,
        atTime=None,
        waitMinutes=0,
        waitHours=0,
        queue=None,
        memory=90,
        memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=True,
        explicitPrecision=DOUBLE,
        nodalOutputPrecision=FULL,
        echoPrint=OFF,
        modelPrint=OFF,
        contactPrint=OFF,
        historyPrint=OFF,
        userSubroutine='',
        scratch='',
        resultsFormat=ODB,
        parallelizationMethodExplicit=DOMAIN,
        numDomains=1,
        activateLoadBalancing=False,
        multiprocessingMode=DEFAULT,
        numCpus=1,
        numGPUs=0
    )
    mdb.jobs[job_name].writeInput(consistencyChecking=OFF)


def stp_to_inp_pipeline(stp_path, output_inp_path,
                        array_num=(1, 1, 1),
                        is_8th=True,
                        mesh_size=0.4):
    """Integrated: STP → array + mesh → export inp, return inp path"""
    if isinstance(array_num, (list, tuple)):
        if len(array_num) != 3:
            raise ValueError("array_num must be length-3 [nx,ny,nz]")
        arrayX_num, arrayY_num, arrayZ_num = map(int, array_num)
    else:
        n = int(array_num)
        arrayX_num = arrayY_num = arrayZ_num = n

    work_path = os.path.dirname(os.path.abspath(output_inp_path))
    if not os.path.exists(work_path):
        os.makedirs(work_path)
    os.chdir(work_path)

    model_name = 'Model-1'
    mdb.Model(name=model_name)

    if is_8th:
        part_name = 'Part-8th'
        x_max, y_max, z_max, x_mid, y_mid, z_mid = import_stp_to_part(
            stp_path, model_name, part_name
        )
        x_center = 2 * x_mid - x_max
        y_center = 2 * y_mid - y_max
        z_center = 2 * z_mid - z_max
        x_len = 2 * (x_max - x_center)
        y_len = 2 * (y_max - y_center)
        z_len = 2 * (z_max - z_center)
    else:
        part_name = 'originStructure'
        x_max, y_max, z_max, x_mid, y_mid, z_mid = import_stp_to_part(
            stp_path, model_name, part_name
        )
        x_center, y_center, z_center = x_mid, y_mid, z_mid
        x_len = 2 * (x_max - x_mid)
        y_len = 2 * (y_max - y_mid)
        z_len = 2 * (z_max - z_mid)

        build_box_for_eighth(model_name, x_max, x_center, y_max, y_center, z_max, z_center)
        cut_to_eighth(model_name, part_name)

    mesh_eighth(model_name, part_name='Part-8th', seed_size=mesh_size)
    mirror_full_cell(model_name, x_center, y_center, z_center)
    build_array(model_name, arrayX_num, arrayY_num, arrayZ_num,
                x_len, y_len, z_len,
                x_center, y_center, z_center)
    clean_intermediate_instances(model_name)

    job_name = os.path.splitext(os.path.basename(output_inp_path))[0]
    write_inp_from_model(model_name, job_name)
    print("INP written to:", output_inp_path)
    return output_inp_path


# ============================================================
# 2. INP → steady-state heat transfer analysis → ODB
# ============================================================

def import_inp_to_model(inp_path, model_name=None):
    """Import model from inp and return model_name"""
    filename_with_ext = os.path.basename(inp_path)
    default_name = os.path.splitext(filename_with_ext)[0]
    if model_name is None:
        model_name = default_name
    mdb.ModelFromInputFile(name=model_name, inputFileName=inp_path)
    return model_name


def change_mesh_to_heat(model_name, part_name,
                        elem_code=DC3D10, elem_library=STANDARD):
    """Change the element type of a given Part to heat-transfer elements"""
    p = mdb.models[model_name].parts[part_name]
    elem_type = mesh.ElemType(elemCode=elem_code, elemLibrary=elem_library)
    all_elems = p.elements
    p.setElementType(regions=(all_elems,), elemTypes=(elem_type,))


def set_boundary_sets(model_name, instance_name, mesh_range=1e-4):
    """Automatically create XTOP/XBOT/YTOP/.../ALLBOUNDARY node sets based on bounding box"""
    a = mdb.models[model_name].rootAssembly
    inst = a.instances[instance_name]
    nodes = inst.nodes

    x_vals = [n.coordinates[0] for n in nodes]
    y_vals = [n.coordinates[1] for n in nodes]
    z_vals = [n.coordinates[2] for n in nodes]

    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)
    z_min, z_max = min(z_vals), max(z_vals)

    # Use UPPERCASE set names, to match cfg: set_name = 'YTOP'
    boundary_sets = ['XTOP', 'XBOT', 'YTOP', 'YBOT', 'ZTOP', 'ZBOT', 'ALLBOUNDARY']

    for bname in boundary_sets:
        if bname.upper() == 'ALLBOUNDARY':
            # whole boundary
            nodes_on_box = nodes
        else:
            # init bounding box dict with lowercase axis keys
            coor = {
                'xMin': x_min, 'xMax': x_max,
                'yMin': y_min, 'yMax': y_max,
                'zMin': z_min, 'zMax': z_max
            }

            # axis: 'x'/'y'/'z'
            axis = bname[0].lower()
            name_upper = bname.upper()

            if 'TOP' in name_upper:
                coor['%sMin' % axis] = coor['%sMax' % axis] - mesh_range
                coor['%sMax' % axis] = coor['%sMax' % axis] + mesh_range
            elif 'BOT' in name_upper:
                coor['%sMax' % axis] = coor['%sMin' % axis] + mesh_range
                coor['%sMin' % axis] = coor['%sMin' % axis] - mesh_range

            nodes_on_box = nodes.getByBoundingBox(**coor)

        print("Get %s: %d nodes" % (bname, len(nodes_on_box)))
        a.Set(nodes=nodes_on_box, name=bname)


def apply_temperature_bcs(model_name, direction,
                          Thot=30.0, Tcold=25.0, Tinit=25.0):
    """
    Apply temperature BCs on direction TOP/BOT and initial temperature on ALLBOUNDARY.
    direction: 'x' / 'y' / 'z'
    """
    a = mdb.models[model_name].rootAssembly
    a.regenerate()

    # create step
    mdb.models[model_name].HeatTransferStep(
        name='Step-1',
        previous='Initial',
        response=STEADY_STATE,
        amplitude=RAMP
    )

    # set names are uppercase: XTOP, XBOT, ...
    dir_upper = direction.upper()   # 'y' -> 'Y'

    # TOP BC
    region_top = a.sets[dir_upper + 'TOP']
    mdb.models[model_name].TemperatureBC(
        name='BC-' + dir_upper + 'TOP',
        createStepName='Step-1',
        region=region_top,
        fixed=OFF,
        distributionType=UNIFORM,
        fieldName='',
        magnitude=Thot,
        amplitude=UNSET
    )

    # BOTTOM BC
    region_bot = a.sets[dir_upper + 'BOT']
    mdb.models[model_name].TemperatureBC(
        name='BC-' + dir_upper + 'BOT',
        createStepName='Step-1',
        region=region_bot,
        fixed=OFF,
        distributionType=UNIFORM,
        fieldName='',
        magnitude=Tcold,
        amplitude=UNSET
    )

    # Initial temperature on ALLBOUNDARY
    region_init = a.sets['ALLBOUNDARY']
    mdb.models[model_name].Temperature(
        name='Predefined Field-1',
        createStepName='Initial',
        region=region_init,
        distributionType=UNIFORM,
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
        magnitudes=Tinit,
    )



def assign_material(model_name, part_name,
                    k=180.0, rho=2.7e-9, cp=9.0e8,
                    mat_name='Material-1', sec_name='Section-1'):
    """Define material & section, and assign to all elements (units follow original script)"""
    model = mdb.models[model_name]
    model.Material(name=mat_name)
    model.materials[mat_name].Conductivity(table=((k,),))
    model.materials[mat_name].Density(table=((rho,),))
    model.materials[mat_name].SpecificHeat(table=((cp,),))

    model.HomogeneousSolidSection(name=sec_name, material=mat_name, thickness=None)

    p = model.parts[part_name]
    elems = p.elements
    region = p.Set(elements=elems, name='Set-AllElems')
    p.SectionAssignment(
        region=region,
        sectionName=sec_name,
        offset=0.0,
        offsetType=MIDDLE_SURFACE,
        offsetField='',
        thicknessAssignment=FROM_SECTION
    )


def submit_heat_job(model_name, ncpus=1, queue=None):
    """Create and submit a steady-state heat-transfer Job, return ABSOLUTE odb path"""
    assert ncpus > 0

    job_name = model_name

    # 1) define a local scratch folder under current working directory
    work_dir = os.getcwd()
    scratch_dir = os.path.join(work_dir, "abaqus_scratch")
    if not os.path.exists(scratch_dir):
        os.makedirs(scratch_dir)

    print("Current working directory for job:", work_dir)
    print("Scratch directory for job:", scratch_dir)
    print("%s  Solving model '%s' with %d CPUs..." %
          (time.ctime(), model_name, ncpus))

    if job_name in mdb.jobs.keys():
        del mdb.jobs[job_name]

    mdb.Job(
        name=job_name,
        model=model_name,
        description='steady-state heat transfer',
        type=ANALYSIS,
        atTime=None,
        waitMinutes=0,
        waitHours=0,
        queue=queue,
        memory=90,
        memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=True,
        explicitPrecision=DOUBLE,
        nodalOutputPrecision=FULL,
        echoPrint=OFF,
        modelPrint=OFF,
        contactPrint=OFF,
        historyPrint=OFF,
        userSubroutine='',
        scratch=scratch_dir,         
        resultsFormat=ODB,
        multiprocessingMode=DEFAULT,
        numCpus=ncpus,
        numDomains=ncpus,
        numGPUs=0
    )
    mdb.jobs[job_name].submit(consistencyChecking=OFF)
    mdb.jobs[job_name].waitForCompletion()

    status = mdb.jobs[job_name].status
    print("%s  Job '%s' finished with status: %s" %
          (time.ctime(), job_name, status))

    # if status != COMPLETED:
    #     raise RuntimeError(
    #         "Abaqus job '%s' did NOT complete successfully (status=%s). "
    #         "Check %s.dat/.log/.msg for details."
    #         % (job_name, status, job_name)
    #     )

    odb_path = os.path.abspath(job_name + '.odb')
    return odb_path


def inp_to_odb_pipeline(inp_path,
                        direction='y',
                        part_name='FINALINSTANCE',
                        instance_name='FINALINSTANCE-1',
                        mesh_range=1e-4,
                        Thot=30.0,
                        Tcold=25.0,
                        Tinit=25.0,
                        k=180.0,
                        rho=2.7e-9,
                        cp=9.0e8,
                        ncpus=1,
                        queue=None):
    """Integrated: INP → heat-transfer model → solve → return odb path"""
    
    filename_with_ext = os.path.basename(inp_path)
    default_name = os.path.splitext(filename_with_ext)[0]
    model_name = default_name + '_' + direction
    
    import_inp_to_model(inp_path, model_name)
    print("Imported model:", model_name)
    
    change_mesh_to_heat(model_name, part_name)
    set_boundary_sets(model_name, instance_name, mesh_range=mesh_range)
    apply_temperature_bcs(model_name, direction, Thot=Thot, Tcold=Tcold, Tinit=Tinit)
    assign_material(model_name, part_name, k=k, rho=rho, cp=cp)
    odb_path = submit_heat_job(model_name, ncpus=ncpus, queue=queue)
    return odb_path


# ============================================================
# 3. ODB → average HFL → keff
# ============================================================

def get_average_HFL_gui(odb_path, set_name, direction):
    print('first_checkpoint')
    """Extract the average HFL component on a given node set from odb using GUI session (unweighted average)"""
    
    value_map = {'x': 'HFL1', 'y': 'HFL2', 'z': 'HFL3'}
    HFL_name = value_map.get(direction)
    if HFL_name is None:
        raise ValueError("direction must be 'x','y','z'")
    
    o1 = session.openOdb(name=odb_path)
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    
    odb = session.odbs[odb_path]
    
    session.xyDataListFromField(
        odb=odb, 
        outputPosition=NODAL, 
        variable=(('HFL', INTEGRATION_POINT, ((COMPONENT, HFL_name), )), ), 
        nodeSets=(set_name, )
    )
    
    xy1 = session.xyDataObjects
    node_names = xy1.keys()
    node_num = len(node_names)
    
    HFL_sum = 0
    for node_name in node_names:
        HFL_node = xy1[node_name].data[-1][-1]
        HFL_sum += HFL_node
        del xy1[node_name]
    
    avg_hfl = HFL_sum / node_num if node_num > 0 else 0.0
    
    print("average HFL (gui, simple):", avg_hfl)
    print("point number:", node_num)
    
    return float(avg_hfl)


def compute_keff_from_HFL(avg_HFL, L, deltaT):
    """Compute keff = q * L / ΔT from average heat flux avg_HFL, length L and temperature difference deltaT"""
    if deltaT == 0:
        raise ZeroDivisionError("deltaT cannot be zero")
    return avg_HFL * L / deltaT


# ============================================================
# 4. Clean intermediate files
# ============================================================

def cleanup_job_files(job_basename, keep_odb=True):
    """
    Delete intermediate files generated by Abaqus Job, optionally keep the odb file.
    job_basename: job name without extension (usually = model_name / inp basename)
    """
    exts = [
        '.com', '.dat', '.ipm', '.log', '.msg',
        '.prt', '.res', '.sta', '.mdl', '.stt',
        '.sim', '.lck'
    ]
    for ext in exts:
        f = job_basename + ext
        if os.path.exists(f):
            try:
                os.remove(f)
                print("Deleted:", f)
            except Exception as e:
                print("Cannot delete %s: %s" % (f, e))

    if not keep_odb:
        odb_file = job_basename + '.odb'
        if os.path.exists(odb_file):
            try:
                os.remove(odb_file)
                print("Deleted:", odb_file)
            except Exception as e:
                print("Cannot delete %s: %s" % (odb_file, e))


# ============================================================
# 5. New function: Run simulations in all three directions and collect results
# ============================================================

def run_all_directions(stp_path,
                       output_inp_path,
                       array_num=(1, 1, 1),
                       is_8th=True,
                       mesh_size=0.4,
                       part_name='FINALINSTANCE',
                       instance_name='FINALINSTANCE-1',
                       mesh_range=1e-4,
                       Thot=30.0,
                       Tcold=25.0,
                       Tinit=25.0,
                       k=180.0,
                       rho=2.7e-9,
                       cp=9.0e8,
                       L=10.0,
                       deltaT=5.0,
                       ncpus=1,
                       queue=None,
                       keep_odb=True,
                       result_csv='keff_results_all_directions.csv'):
    """
    Modification: Run simulations in three directions, calculate and output 
    heat flux density and thermal conductivity in x, y, z directions
    
    Parameters:
    -----------
    result_csv : str
        Output CSV file path, will contain HFL1, HFL2, HFL3 and corresponding thermal conductivities
    """
    
    # First convert STP to INP (do this only once)
    inp_path = stp_to_inp_pipeline(
        stp_path=stp_path,
        output_inp_path=output_inp_path,
        array_num=array_num,
        is_8th=is_8th,
        mesh_size=mesh_size
    )
    
    # Initialize results dictionary
    results = {}
    
    # Modification: Loop through simulations in three directions
    directions = ['x', 'y', 'z']
    
    for direction in directions:
        print("\n" + "="*60)
        print("Starting simulation in {} direction".format(direction))
        print("="*60)
        
        # Modification: Set corresponding boundary condition sets
        if direction == 'x':
            set_name = 'XTOP'
        elif direction == 'y':
            set_name = 'YTOP'
        else:  # direction == 'z'
            set_name = 'ZTOP'
        
        # Run thermal analysis in this direction
        odb_path = inp_to_odb_pipeline(
            inp_path=inp_path,
            direction=direction,
            part_name=part_name,
            instance_name=instance_name,
            mesh_range=mesh_range,
            Thot=Thot,
            Tcold=Tcold,
            Tinit=Tinit,
            k=k,
            rho=rho,
            cp=cp,
            ncpus=ncpus,
            queue=queue
        )
        
        # Extract heat flux density in this direction
        avg_HFL = get_average_HFL_gui(
            odb_path=odb_path,
            set_name=set_name,
            direction=direction
        )
        
        # Calculate thermal conductivity in this direction
        keff = compute_keff_from_HFL(avg_HFL, L=L, deltaT=deltaT)
        
        # Store results
        results[direction] = {
            'avg_HFL': avg_HFL,
            'keff': keff,
            'odb_path': odb_path
        }
        
        print("{} direction result: avg_HFL={:.6e}, keff={:.6f}".format(direction, avg_HFL, keff))
        
        # Clean up files (optionally keep ODB)
        job_base = os.path.splitext(os.path.basename(inp_path))[0] + '_' + direction
        cleanup_job_files(job_base, keep_odb=keep_odb)
    
    # Output results to CSV file
    # Modification: Use new column names HFL1, HFL2, HFL3 and corresponding thermal conductivities
    try:
        header = "stp_file,HFL1_x,keff1_x,HFL2_y,keff2_y,HFL3_z,keff3_z,L\n"
        line = "%s,%g,%g,%g,%g,%g,%g,%g\n" % (
            os.path.basename(stp_path),
            results['x']['avg_HFL'], results['x']['keff'],
            results['y']['avg_HFL'], results['y']['keff'],
            results['z']['avg_HFL'], results['z']['keff'],
            L
        )
        
        # Modification: Use new CSV filename
        result_csv = '/public/home/nieqi01/hzz/abaqus/test_single_kave/keff_results_all_directions.csv'
        
        if (not os.path.exists(result_csv)) or os.path.getsize(result_csv) == 0:
            with open(result_csv, 'w') as f:
                f.write(header)
                f.write(line)
        else:
            with open(result_csv, 'a') as f:
                f.write(line)
        
        print("\n" + "="*60)
        print("All direction simulations completed!")
        print("Results saved to: {}".format(result_csv))
        print("="*60)
        
        # Print summary results
        print("\nSummary results:")
        print("x direction: HFL1 = {:.6e}, keff = {:.6f}".format(results['x']['avg_HFL'], results['x']['keff']))
        print("y direction: HFL2 = {:.6e}, keff = {:.6f}".format(results['y']['avg_HFL'], results['y']['keff']))
        print("z direction: HFL3 = {:.6e}, keff = {:.6f}".format(results['z']['avg_HFL'], results['z']['keff']))
        
    except Exception as e:
        print("Error writing to CSV file:", e)
    
    return results


# ============================================================
# 6. Configuration section: Control your task here
# ============================================================

if __name__ == '__main__':
    # Allow overriding via environment variables (for Slurm on cluster)
    default_stp = "/public/home/nieqi01/hzz/abaqus/test_single_kave/137777.stp"
    default_inp = "/public/home/nieqi01/hzz/abaqus/test_single_kave/137777-symmetryMesh.inp"

    stp_path = os.getenv("BPMESH_INPUT_MODEL", default_stp)
    output_inp_path = os.getenv("BPMESH_OUTPUT_MESH", default_inp)
    is_8th = bool(int(os.getenv("BPMESH_IS_8TH", '0')))  # 1 = input is a 1/8 unit cell

    # Array size, can be changed to [2,2,2], etc.
    array_num = [1, 1, 1]

    # Modification: No need to specify a single direction, as all directions will be calculated
    L = 10.0   # e.g., 10 mm
    Thot = 30.0
    Tcold = 25.0
    deltaT = Thot - Tcold

    # Default naming for BCC 137777 case: part / instance
    part_name = 'FINALINSTANCE'
    instance_name = 'FINALINSTANCE-1'

    # Resources and outputs
    ncpus = 1
    queue = None
    keep_odb = True
    
    # Modification: Call the new function to run all directions
    results = run_all_directions(
        stp_path=stp_path,
        output_inp_path=output_inp_path,
        array_num=array_num,
        is_8th=is_8th,
        mesh_size=0.4,
        part_name=part_name,
        instance_name=instance_name,
        mesh_range=1e-4,
        Thot=Thot,
        Tcold=Tcold,
        Tinit=Tcold,   # initial temperature set equal to cold side
        k=180.0,
        rho=2.7e-9,
        cp=9.0e8,
        L=L,
        deltaT=deltaT,
        ncpus=ncpus,
        queue=queue,
        keep_odb=keep_odb
    )