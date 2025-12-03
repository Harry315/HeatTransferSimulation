# -*- coding: utf-8 -*-
"""
Full pipeline: STP → mesh/inp → heat analysis → HFL → keff

步骤：
1. 从 STP 生成 1/8 单胞几何 → 镜像成完整单胞 → 阵列 → 网格 → 导出 inp
2. 从 inp 导入模型，转换为热传导模型，施加边界 & 材料，求解
3. 打开 odb，读取指定节点集上的 HFL，计算平均热流 q
4. 根据 q、L、ΔT 计算有效导热率 keff
5. 删除中间文件（保留 odb 和结果 csv）

说明：
- 适合在 abaqus python nogui 模式下跑（HPC / 曙光场景）
- 若要用 GUI 模式的 HFL 提取，可额外加上 session 相关函数
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
# 1. STP → 1/8 单胞 → 镜像 → 阵列 → 网格 → 导出 inp
# ============================================================

def import_stp_to_part(stp_path, model_name, part_name):
    """读取 STP，创建 Part，并返回最大坐标和中心坐标"""
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
    """画 box（1/8 空间），用于从完整结构切出 1/8 单胞"""
    m = mdb.models[model_name]
    s = m.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(x_center, y_center), point2=(x_max, y_max))

    p = m.Part(name='box', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s, depth=(z_max - z_center))

    s.unsetPrimaryObject()
    del m.sketches['__profile__']


def cut_to_eighth(model_name, part_name):
    """利用 box 和原始结构做布尔运算，得到 Part-8th"""
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
    """对 Part-8th 网格划分"""
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
    """1/8 → 1/4 → 1/2 → 整单胞 mirror-finished"""
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

    # 整体 + XY mirror
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
    """用 mirror-finished 做三向阵列，合并为 finalInstance"""
    m = mdb.models[model_name]
    a = m.rootAssembly

    # XY 阵列
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

    # Z 阵列
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

    # 平移居中
    a.translate(
        instanceList=('finalInstance-1',),
        vector=(
            x_center - arrayX_num * x_len * 0.5,
            y_center - arrayY_num * y_len * 0.5,
            z_center - arrayZ_num * z_len * 0.5
        )
    )


def clean_intermediate_instances(model_name):
    """删除中间实例，只保留 finalInstance-1"""
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
    """建 job 并写 inp（不求解）"""
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
        explicitPrecision=SINGLE,
        nodalOutputPrecision=SINGLE,
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
    """整合：STP → 阵列网格 → 导出 inp，返回 inp 路径"""
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
# 2. INP → 热传导分析（稳态）→ ODB
# ============================================================

def import_inp_to_model(inp_path, model_name=None):
    """从 inp 导入模型，返回 model_name"""
    filename_with_ext = os.path.basename(inp_path)
    default_name = os.path.splitext(filename_with_ext)[0]
    if model_name is None:
        model_name = default_name
    mdb.ModelFromInputFile(name=model_name, inputFileName=inp_path)
    return model_name


def change_mesh_to_heat(model_name, part_name,
                        elem_code=DC3D10, elem_library=STANDARD):
    """把指定 part 的单元类型改成热传导单元"""
    p = mdb.models[model_name].parts[part_name]
    elem_type = mesh.ElemType(elemCode=elem_code, elemLibrary=elem_library)
    all_elems = p.elements
    p.setElementType(regions=(all_elems,), elemTypes=(elem_type,))


def set_boundary_sets(model_name, instance_name, mesh_range=1e-4):
    """根据包络盒自动生成 xtop/xbot/ytop/.../allBoundary 节点集"""
    a = mdb.models[model_name].rootAssembly
    inst = a.instances[instance_name]
    nodes = inst.nodes

    x_vals = [n.coordinates[0] for n in nodes]
    y_vals = [n.coordinates[1] for n in nodes]
    z_vals = [n.coordinates[2] for n in nodes]

    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)
    z_min, z_max = min(z_vals), max(z_vals)

    boundary_sets = ['xtop', 'xbot', 'ytop', 'ybot', 'ztop', 'zbot', 'allBoundary']

    for bname in boundary_sets:
        if bname == 'allBoundary':
            nodes_on_box = nodes
        else:
            coor = {
                'xMin': x_min, 'xMax': x_max,
                'yMin': y_min, 'yMax': y_max,
                'zMin': z_min, 'zMax': z_max
            }
            direction = bname[0]  # x/y/z
            if 'top' in bname:
                coor['%sMin' % direction] = coor['%sMax' % direction] - mesh_range
                coor['%sMax' % direction] = coor['%sMax' % direction] + mesh_range
            elif 'bot' in bname:
                coor['%sMax' % direction] = coor['%sMin' % direction] + mesh_range
                coor['%sMin' % direction] = coor['%sMin' % direction] - mesh_range

            nodes_on_box = nodes.getByBoundingBox(**coor)

        print("Get %s: %d nodes" % (bname, len(nodes_on_box)))
        a.Set(nodes=nodes_on_box, name=bname)


def apply_temperature_bcs(model_name, direction,
                          Thot=30.0, Tcold=25.0, Tinit=25.0):
    """在 direction(top/bot) 上施温度边界，并在 allBoundary 上施初始温度"""
    a = mdb.models[model_name].rootAssembly
    a.regenerate()

    mdb.models[model_name].HeatTransferStep(
        name='Step-1',
        previous='Initial',
        response=STEADY_STATE,
        amplitude=RAMP
    )

    region_top = a.sets[direction + 'top']
    mdb.models[model_name].TemperatureBC(
        name='BC-' + direction + 'top',
        createStepName='Step-1',
        region=region_top,
        fixed=OFF,
        distributionType=UNIFORM,
        fieldName='',
        magnitude=Thot,
        amplitude=UNSET
    )

    region_bot = a.sets[direction + 'bot']
    mdb.models[model_name].TemperatureBC(
        name='BC-' + direction + 'bot',
        createStepName='Step-1',
        region=region_bot,
        fixed=OFF,
        distributionType=UNIFORM,
        fieldName='',
        magnitude=Tcold,
        amplitude=UNSET
    )

    region_init = a.sets['allBoundary']
    mdb.models[model_name].Temperature(
        name='Predefined Field-1',
        createStepName='Initial',
        region=region_init,
        distributionType=UNIFORM,
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
        magnitudes=(Tinit,)
    )


def assign_material(model_name, part_name,
                    k=180.0, rho=2.7e-9, cp=9.0e8,
                    mat_name='Material-1', sec_name='Section-1'):
    """定义材料 & 截面，并赋给所有单元（单位制沿用原脚本）"""
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
    """创建并提交稳态热传导 Job，返回 odb 路径"""
    assert ncpus > 0
    print("%s  Solving model '%s' with %d CPUs..." %
          (time.ctime(), model_name, ncpus))

    job_name = model_name
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
        explicitPrecision=SINGLE,
        nodalOutputPrecision=SINGLE,
        echoPrint=OFF,
        modelPrint=OFF,
        contactPrint=OFF,
        historyPrint=OFF,
        userSubroutine='',
        scratch='',
        resultsFormat=ODB,
        multiprocessingMode=DEFAULT,
        numCpus=ncpus,
        numDomains=ncpus,
        numGPUs=0
    )
    mdb.jobs[job_name].submit(consistencyChecking=OFF)
    mdb.jobs[job_name].waitForCompletion()
    print("%s  Job '%s' finished." % (time.ctime(), job_name))

    odb_path = job_name + '.odb'
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
    """整合：INP → 热模型 → 求解 → 返回 odb 路径"""
    model_name = import_inp_to_model(inp_path)
    print("Imported model:", model_name)

    change_mesh_to_heat(model_name, part_name)
    set_boundary_sets(model_name, instance_name, mesh_range=mesh_range)
    apply_temperature_bcs(model_name, direction, Thot=Thot, Tcold=Tcold, Tinit=Tinit)
    assign_material(model_name, part_name, k=k, rho=rho, cp=cp)
    odb_path = submit_heat_job(model_name, ncpus=ncpus, queue=queue)
    return odb_path


# ============================================================
# 3. ODB → 平均 HFL → keff
# ============================================================

def get_average_HFL_nogui(odb_path, set_name, direction):
    """从 odb 中提取指定节点集上的 HFL 分量平均值（不加权）"""
    value_map = {'x': 0, 'y': 1, 'z': 2}
    comp_idx = value_map.get(direction)
    if comp_idx is None:
        raise ValueError("direction must be 'x','y','z'")

    odb = openOdb(path=odb_path, readOnly=True)

    # 最后一步最后一帧
    last_step = odb.steps[odb.steps.keys()[-1]]
    last_frame = last_step.frames[-1]

    node_set = odb.rootAssembly.nodeSets[set_name]
    node_labels = [n.label for n in node_set.nodes[0]]

    hfl_nodal = last_frame.fieldOutputs['HFL'].getSubset(position=ELEMENT_NODAL)
    hfl_region = hfl_nodal.getSubset(region=node_set)

    last_val = {}
    for v in hfl_region.values:
        last_val[v.nodeLabel] = v.data[comp_idx]

    hfl_array = np.array([last_val[lbl] for lbl in node_labels if lbl in last_val])
    avg_hfl = float(hfl_array.mean())
    odb.close()

    print("average HFL (nogui, simple):", avg_hfl)
    print("point number:", len(hfl_array))
    return avg_hfl


def compute_keff_from_HFL(avg_HFL, L, deltaT):
    """根据平均热流 avg_HFL、长度 L、ΔT 计算有效导热率 keff = q * L / ΔT"""
    if deltaT == 0:
        raise ZeroDivisionError("deltaT cannot be zero")
    return avg_HFL * L / deltaT


# ============================================================
# 4. 清理中间文件
# ============================================================

def cleanup_job_files(job_basename, keep_odb=True):
    """
    删除 Abaqus Job 生成的中间文件，只保留 odb（可选）
    job_basename: 不带后缀的 job 名（通常 = model_name / inp basename）
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
# 5. 总控制：一条龙 pipeline
# ============================================================

def full_pipeline(stp_path,
                  output_inp_path,
                  array_num=(1, 1, 1),
                  is_8th=True,
                  mesh_size=0.4,
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
                  L=10.0,
                  deltaT=5.0,
                  ncpus=1,
                  queue=None,
                  set_name='YTOP',
                  keep_odb=True,
                  result_csv='keff_results.csv'):
    """
    一站式：STP → INP → ODB → HFL → keff
    - L: 试样在传热方向上的长度（单位与模型一致）
    - deltaT: top 与 bot 的温差 = Thot - Tcold
    """

    # 1. STP → INP
    inp_path = stp_to_inp_pipeline(
        stp_path=stp_path,
        output_inp_path=output_inp_path,
        array_num=array_num,
        is_8th=is_8th,
        mesh_size=mesh_size
    )

    # 2. INP → ODB（稳态热传导）
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

    # 3. ODB → 平均 HFL → keff
    avg_HFL = get_average_HFL_nogui(
        odb_path=odb_path,
        set_name=set_name,
        direction=direction
    )
    keff = compute_keff_from_HFL(avg_HFL, L=L, deltaT=deltaT)
    print("Effective thermal conductivity keff =", keff)

    # 4. 写结果到 csv（追加）
    try:
        line = "%s,%s,%g,%g,%g\n" % (
            os.path.basename(stp_path),
            os.path.basename(odb_path),
            avg_HFL,
            keff,
            L
        )
        header = "stp_file,odb_file,avg_HFL,keff,L\n"
        if (not os.path.exists(result_csv)) or os.path.getsize(result_csv) == 0:
            with open(result_csv, 'w') as f:
                f.write(header)
                f.write(line)
        else:
            with open(result_csv, 'a') as f:
                f.write(line)
        print("Result appended to:", result_csv)
    except Exception as e:
        print("Cannot write csv:", e)

    # 5. 清理中间文件
    job_base = os.path.splitext(os.path.basename(output_inp_path))[0]
    cleanup_job_files(job_base, keep_odb=keep_odb)

    return keff, avg_HFL, odb_path


# ============================================================
# 6. 配置区：只改这里来控制你的任务
# ============================================================

if __name__ == '__main__':
    # 允许用环境变量覆盖，方便在曙光上用 Slurm 传参
    default_stp = r'C:\Users\Administrator\Desktop\20251022\calHFL\137777.stp'
    default_inp = r'C:\Users\Administrator\Desktop\20251022\calHFL\137777-symmetryMesh.inp'

    stp_path = os.getenv("BPMESH_INPUT_MODEL", default_stp)
    output_inp_path = os.getenv("BPMESH_OUTPUT_MESH", default_inp)
    is_8th = bool(int(os.getenv("BPMESH_IS_8TH", '1')))  # 1 = 输入是 1/8 单胞

    # 阵列数量，可以改成 [2,2,2] 之类
    array_num = [1, 1, 1]

    # 传热方向 & 试样长度（L 建议和 BCC 单胞尺寸 * 阵列数一致）
    direction = 'y'
    L = 10.0   # 例如 10 mm
    Thot = 30.0
    Tcold = 25.0
    deltaT = Thot - Tcold

    # BCC 137777 那套默认：part / instance 命名
    part_name = 'FINALINSTANCE'
    instance_name = 'FINALINSTANCE-1'
    set_name = 'YTOP'  # 对应 direction='y' 顶面

    # 资源与输出
    ncpus = 1
    queue = None
    keep_odb = True
    result_csv = 'keff_results.csv'

    full_pipeline(
        stp_path=stp_path,
        output_inp_path=output_inp_path,
        array_num=array_num,
        is_8th=is_8th,
        mesh_size=0.4,
        direction=direction,
        part_name=part_name,
        instance_name=instance_name,
        mesh_range=1e-4,
        Thot=Thot,
        Tcold=Tcold,
        Tinit=Tcold,   # 初始温度设成冷端温度
        k=180.0,
        rho=2.7e-9,
        cp=9.0e8,
        L=L,
        deltaT=deltaT,
        ncpus=ncpus,
        queue=queue,
        set_name=set_name,
        keep_odb=keep_odb,
        result_csv=result_csv
    )
