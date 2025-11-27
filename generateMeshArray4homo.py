from abaqus import *
from abaqusConstants import *
import __main__
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
import time
import multiprocessing
import os
import ctypes
from caeModules import *
from driverUtils import executeOnCaeStartup
import os
import csv


def main(stp_path,output_name,array_num,is_8th=True):
    # if array_num is list or tuple, it can only have 3 elements, which indicate the array number in xyz direction
    # if array_num is number, it mean the array number is the same in xyz direction as array_num
    
    if isinstance(array_num,(list,tuple)):
        if len(array_num) != 3:                
            raise ValueError("Elements number in array_num must be three")
        arrayX_num,arrayY_num,arrayZ_num=array_num[0],array_num[1],array_num[2]
    else:
        arrayX_num,arrayY_num,arrayZ_num=int(array_num),int(array_num),int(array_num)



    print("inputPath:",stp_path)

    work_path=os.path.dirname(os.path.abspath(output_name))
    # check if folder exists\
    if not os.path.exists(work_path):
        os.makedirs(work_path)
    
    # full_path = os.path.join(job_path, folder_path)
    os.chdir(work_path)
    path = os.getcwd()
    print("Current working directory:", path)


    # step = mdb.openStep(os.path.join(model_dir, files[i]), scaleFromFile=OFF)
    # # mdb.models['Model-1'].PartFromGeometryFile(name='e1', geometryFile=step, 
    # #     combine=True, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    # mdb.models['Model-1'].PartFromGeometryFile(name='e1', geometryFile=step, 
    #     combine=True, retainBoundary=True, mergeSolidRegions=True, 
    #     dimensionality=THREE_D, type=DEFORMABLE_BODY)    
    # p = mdb.models['Model-1'].parts['e1']

    # x_coords = [vertex.pointOn[0][0] for vertex in p.vertices]

    # x_min = round(min(x_coords), 3) 
    # x_max =round(max(x_coords), 3) 
    # x_middle=(x_max+x_min)/2



    # print('min:',x_min)
    # print(x_max)        


    
    def importstp(stp_path,part_name):

        step = mdb.openStep((stp_path), scaleFromFile=OFF)
        mdb.models['Model-1'].PartFromGeometryFile(name=part_name, geometryFile=step, 
            combine=True, mergeSolidRegions=True, 
            dimensionality=THREE_D, type=DEFORMABLE_BODY)  


        p = mdb.models['Model-1'].parts[part_name]

        x_coords = [vertex.pointOn[0][0] for vertex in p.vertices]

        x_min = round(min(x_coords), 3) 
        x_max =round(max(x_coords), 3) 
        x_middle=(x_max+x_min)/2

        y_coords = [vertex.pointOn[0][1] for vertex in p.vertices]

        y_min = round(min(y_coords), 3) 
        y_max =round(max(y_coords), 3) 
        y_middle=(y_max+y_min)/2

        z_coords = [vertex.pointOn[0][2] for vertex in p.vertices]

        z_min = round(min(z_coords), 3) 
        z_max =round(max(z_coords), 3) 
        z_middle=(z_max+z_min)/2
        return x_max,y_max,z_max,x_middle,y_middle,z_middle


    def buildbox(x_max,x_middle,y_max,y_middle,z_max,z_middle):


        # model_name = 'Model-done'
        # if model_name in mdb.models.keys():
        #     del mdb.models[model_name]

        s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
        #g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
        s.setPrimaryObject(option=STANDALONE)
        s.rectangle(point1=(x_middle, y_middle), point2=(x_max, y_max))
        p = mdb.models['Model-1'].Part(name='box', dimensionality=THREE_D, 
            type=DEFORMABLE_BODY)
        p = mdb.models['Model-1'].parts['box']
        p.BaseSolidExtrude(sketch=s, depth=(z_max-z_middle))
        s.unsetPrimaryObject()
        p = mdb.models['Model-1'].parts['box']
        del mdb.models['Model-1'].sketches['__profile__']

    def cut(part_name):
        a = mdb.models['Model-1'].rootAssembly
        a.DatumCsysByDefault(CARTESIAN)
        p = mdb.models['Model-1'].parts['box']
        a.Instance(name='box-1', part=p, dependent=ON)
        p = mdb.models['Model-1'].parts[part_name]
        instance_name=part_name+'-1'
        a.Instance(name=instance_name, part=p, dependent=ON)
        a = mdb.models['Model-1'].rootAssembly

        a = mdb.models['Model-1'].rootAssembly
        print('geometryValidity',a.instances[instance_name].part.geometryValidity) 

        a.InstanceFromBooleanCut(name='Part-boxSub8th', 
            instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances['box-1'], 
            cuttingInstances=(a.instances[instance_name], ), originalInstances=SUPPRESS)
        
        a = mdb.models['Model-1'].rootAssembly
        p = mdb.models['Model-1'].parts['box']
        a.Instance(name='box-2', part=p, dependent=ON)
        a = mdb.models['Model-1'].rootAssembly
        a.InstanceFromBooleanCut(name='Part-8th', 
            instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances['box-2'], 
            cuttingInstances=(a.instances['Part-boxSub8th-1'], ), 
            originalInstances=SUPPRESS)
        #p = mdb.models['Model-1'].parts['e1']
    
    def wangge():
        
        p = mdb.models['Model-1'].parts['Part-8th']
        c = p.cells
        pickedRegions = c.getSequenceFromMask(mask=('[#1 ]', ), )
        p.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
        elemType1 = mesh.ElemType(elemCode=C3D20R)
        elemType2 = mesh.ElemType(elemCode=C3D15)
        elemType3 = mesh.ElemType(elemCode=C3D10)
        p = mdb.models['Model-1'].parts['Part-8th']
        c = p.cells
        cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        pickedRegions =(cells, )
        p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
            elemType3))
        p = mdb.models['Model-1'].parts['Part-8th']
        p.seedPart(size=0.4, deviationFactor=0.08, minSizeFactor=0.01)
        # elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD)
        # elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
        # elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD, 
            # secondOrderAccuracy=OFF, distortionControl=DEFAULT)
        # p = mdb.models['Model-1'].parts['Part-3']
        # c = p.cells
        # cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        # pickedRegions =(cells, )
        # p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
            # elemType3))
        p = mdb.models['Model-1'].parts['Part-8th']
        p.generateMesh()
                  
    def mirror(x_middle,y_middle,z_middle):
        #part:1 and 2 are 8th structure,3 is 4th, 4 is 2nd, mirror-finished is whole structure


        # mdb.models['Model-1'].rootAssembly.instances['Part-3'].makeIndependent()
        p3 = mdb.models['Model-1'].parts['Part-8th']

        p3M = mdb.models['Model-1'].Part(name='Part-8th-MirrorXY', 
            objectToCopy=mdb.models['Model-1'].parts['Part-8th'], 
            compressFeatureList=ON, mirrorPlane=XYPLANE)
        # a = mdb.models['Model-1'].rootAssembly
        a = mdb.models['Model-1'].rootAssembly
        # p = mdb.models['Model-1'].parts['Part-3']
        a.Instance(name='Part-8th-1', part=p3, dependent=ON)
        # p = mdb.models['Model-1'].parts['Part-3-Copy']
        a.Instance(name='Part-8th-MirrorXY-1', part=p3M, dependent=ON)

        # a = mdb.models['Model-1'].rootAssembly
        a.InstanceFromBooleanMerge(name='Part-8th-copy', instances=(a.instances['Part-8th-1'], 
            a.instances['Part-8th-MirrorXY-1'], ), mergeNodes=BOUNDARY_ONLY, 
            nodeMergingTolerance=1e-06, domain=MESH, originalInstances=SUPPRESS)
        #This step didn't show any effect, Part-8th-copy is still Part-8th

        # a1 = mdb.models['Model-1'].rootAssembly
        # a1.translate(instanceList=('Part-3-2', ), vector=(-10.0, -10.0, 
        #     -10.0))
        


        p1 = mdb.models['Model-1'].parts['Part-8th-copy']
        p1M = mdb.models['Model-1'].Part(name='Part-8th-copy-MirrorYZ', 
            objectToCopy=mdb.models['Model-1'].parts['Part-8th-copy'], 
            compressFeatureList=ON, mirrorPlane=YZPLANE)

        a.Instance(name='Part-8th-copy-1', part=p1, dependent=ON)
        
        # p = mdb.models['Model-1'].parts['Part-1-Copy']
        a.Instance(name='Part-8th-copy-MirrorYZ-1', part=p1M, dependent=ON)
        
        a1 = mdb.models['Model-1'].rootAssembly
        a1.translate(instanceList=('Part-8th-copy-1',), vector=(-x_middle, -y_middle, 
            -z_middle))     
        a1 = mdb.models['Model-1'].rootAssembly
        a1.translate(instanceList=('Part-8th-copy-MirrorYZ-1',), vector=(+x_middle, -y_middle, 
            -z_middle))

        # a = mdb.models['Model-1'].rootAssembly
        a.InstanceFromBooleanMerge(name='Part-4th', instances=(a.instances['Part-8th-copy-1'], 
            a.instances['Part-8th-copy-MirrorYZ-1'], ), mergeNodes=BOUNDARY_ONLY, 
            nodeMergingTolerance=1e-06, domain=MESH, originalInstances=SUPPRESS)
        p2 = mdb.models['Model-1'].parts['Part-4th']
        p2M = mdb.models['Model-1'].Part(name='Part-4th-copy-MirrorXZ', 
            objectToCopy=mdb.models['Model-1'].parts['Part-4th'], 
            compressFeatureList=ON, mirrorPlane=XZPLANE)

        a.Instance(name='Part-4th-1', part=p2, dependent=ON)
        # p = mdb.models['Model-1'].parts['Part-2-Copy']
        a.Instance(name='Part-4th-copy-MirrorXZ-1', part=p2M, dependent=ON)
        # a = mdb.models['Model-1'].rootAssembly
        a.InstanceFromBooleanMerge(name='Part-2nd', instances=(a.instances['Part-4th-1'], 
            a.instances['Part-4th-copy-MirrorXZ-1'], ), mergeNodes=BOUNDARY_ONLY, 
            nodeMergingTolerance=1e-06, domain=MESH, originalInstances=SUPPRESS)
        # p = mdb.models['Model-1'].parts['Part-2-Copy']
        # p = mdb.models['Model-1'].parts['Part-3']


        p3 = mdb.models['Model-1'].parts['Part-2nd']
        p3M = mdb.models['Model-1'].Part(name='Part-2nd-copy-MirrorXY-1', 
            objectToCopy=mdb.models['Model-1'].parts['Part-2nd'], 
            compressFeatureList=ON, mirrorPlane=XYPLANE)
        # a = mdb.models['Model-1'].rootAssembly
        a = mdb.models['Model-1'].rootAssembly
        # p = mdb.models['Model-1'].parts['Part-3']
        a.Instance(name='Part-2nd-1', part=p3, dependent=ON)
        # p = mdb.models['Model-1'].parts['Part-3-Copy']
        a.Instance(name='Part-2nd-copy-MirrorXY-1', part=p3M, dependent=ON)

        # a = mdb.models['Model-1'].rootAssembly
        a.InstanceFromBooleanMerge(name='mirror-finished', instances=(a.instances['Part-2nd-1'], 
            a.instances['Part-2nd-copy-MirrorXY-1'], ), mergeNodes=BOUNDARY_ONLY, 
            nodeMergingTolerance=1e-06, domain=MESH, originalInstances=SUPPRESS)
        
    def array(arrayX_num,arrayY_num,arrayZ_num,x_len,y_len,z_len,x_middle,y_middle,z_middle):
        
        a = mdb.models['Model-1'].rootAssembly


        a.LinearInstancePattern(instanceList=('mirror-finished-1', ), direction1=(1.0, 
            0.0, 0.0), direction2=(0.0, 1.0, 0.0), number1=arrayX_num, number2=arrayY_num, 
            spacing1=x_len, spacing2=y_len)
        # a.InstanceFromBooleanMerge(name='array-finished', instances=(a.instances['Part-4-2'], 
        #     a.instances['Part-4-Copy-1'], ), mergeNodes=BOUNDARY_ONLY, 
        #     nodeMergingTolerance=1e-06, domain=MESH, originalInstances=SUPPRESS)

        baseInst  = 'mirror-finished'
        instTuple = tuple(inst for inst in a.instances.values()
                        if inst.name.startswith(baseInst))
        a.InstanceFromBooleanMerge(
            name='array-finishedxy',
            instances=instTuple,          
            mergeNodes=BOUNDARY_ONLY,
            nodeMergingTolerance=1e-06,
            domain=MESH,
            originalInstances=SUPPRESS
        )


        a.LinearInstancePattern(instanceList=('array-finishedxy-1', ), direction1=(1.0, 
            0.0, 0.0), direction2=(0.0, 0.0, 1.0), number1=1, number2=arrayZ_num, 
            spacing1=x_len, spacing2=z_len)

        baseInst  = 'array-finishedxy-1'
        instTuple = tuple(inst for inst in a.instances.values()
                        if inst.name.startswith(baseInst))        
        a.InstanceFromBooleanMerge(
            name='finalInstance',
            instances=instTuple,          
            mergeNodes=BOUNDARY_ONLY,
            nodeMergingTolerance=1e-06,
            domain=MESH,
            originalInstances=SUPPRESS
        )


        a1 = mdb.models['Model-1'].rootAssembly
        print(x_len,y_len,z_len,x_middle,y_middle,z_middle)
        a1.translate(instanceList=('finalInstance-1',), vector=(x_middle-arrayX_num*x_len*0.5, y_middle-arrayY_num*y_len*0.5, 
            z_middle-arrayZ_num*z_len*0.5))    

    def delass():
        a = mdb.models['Model-1'].rootAssembly
        a.regenerate()
        # a.deleteFeatures(('Part-3-1', 'Part-3-2', 'Part-3-Copy-1', 'Part-1-1', 
        #     'Part-1-2', 'Part-1-Copy-1', 'Part-2-1', 'Part-2-2', 'Part-2-Copy-1', ))

        # baseInst  = 'mirror'
        to_del = tuple(inst for inst in a.instances.values()
                        if inst.name.startswith('mirror') or inst.name.lower().startswith('part') or inst.name.lower().startswith('array'))
        for inst in to_del:
            a.deleteFeatures((inst.name,))  
 
    def inpex(job_name):

        if job_name in mdb.jobs.keys():
            del mdb.jobs[job_name]
        else:
            pass
        mdb.Job(name=job_name, model='Model-1', description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, 
            numDomains=1, activateLoadBalancing=False, multiprocessingMode=DEFAULT, 
            numCpus=1, numGPUs=0)
        mdb.jobs[job_name].writeInput(consistencyChecking=OFF)

    modelName = 'Model-1'
    mdb.Model(name=modelName)


    if is_8th:
        part_name='Part-8th'
        x_max,y_max,z_max,x_middle,y_middle,z_middle=importstp(stp_path,part_name)
        x_center=2*x_middle-x_max
        y_center=2*y_middle-y_max
        z_center=2*z_middle-z_max
        x_len=2*(x_max-x_center)
        y_len=2*(y_max-y_center)
        z_len=2*(z_max-z_center)


    else:
        part_name='originStructure'
        x_max,y_max,z_max,x_middle,y_middle,z_middle=importstp(stp_path,part_name)
        x_center=x_middle
        y_center=y_middle
        z_center=z_middle
        x_len=2*(x_max-x_middle)
        y_len=2*(y_max-y_middle)
        z_len=2*(z_max-z_middle)
        buildbox(x_max,x_center,y_max,y_center,z_max,z_center)

        cut(part_name)


    wangge()

    mirror(x_center,y_center,z_center)

    array(arrayX_num,arrayY_num,arrayZ_num,x_len,y_len,z_len,x_center,y_center,z_center)

    delass()
    # maRVE()
    # print('maRVE()')
    outName = os.path.splitext(os.path.basename(output_name))[0]
    inpex(outName)
    print(output_name+' finished!')



stp_path = os.getenv("BPMESH_INPUT_MODEL", r'C:\Users\Administrator\Desktop\20251022\calHFL\137777.stp')   #input stp
output_name = os.getenv("BPMESH_OUTPUT_MESH", r'C:\Users\Administrator\Desktop\20251022\calHFL\137777-symmetryMesh.inp') #output inp
is_8th = int(os.getenv("BPMESH_IS_8TH", '0'))

array_nums=[[1,1,1]]

for number in array_nums:
    main(stp_path,output_name,number,is_8th)


