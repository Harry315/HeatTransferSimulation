# -*- coding: mbcs -*-
# Do not delete the following import lines
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
import os

class stp2mixInstance:
    def __init__(self, stp_path,part_name,mixPartName):
        self.stp_path = stp_path
        self.part_name = part_name
        self.mixPartName=mixPartName

    def importStp(self,stp_path,part_name):
        step = mdb.openStep((stp_path), scaleFromFile=OFF)
        mdb.models['Model-1'].PartFromGeometryFile(name=part_name, geometryFile=step, 
            combine=True, mergeSolidRegions=True, 
            dimensionality=THREE_D, type=DEFORMABLE_BODY)  

        print(part_name)
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
        return x_max,y_max,z_max,x_middle,y_middle,z_middle,x_min,y_min,z_min

    def buildBox(self,x_max,y_max,z_max,x_min,y_min,z_min):
        s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
            sheetSize=200.0)
        s.setPrimaryObject(option=STANDALONE)
        s.rectangle(point1=(x_min,y_min), point2=(x_max,y_max))
        p = mdb.models['Model-1'].Part(name='Part-2', dimensionality=THREE_D, 
            type=DEFORMABLE_BODY)
        p = mdb.models['Model-1'].parts['Part-2']
        p.BaseSolidExtrude(sketch=s, depth=(z_max-z_min))
        s.unsetPrimaryObject()
        p = mdb.models['Model-1'].parts['Part-2']
        del mdb.models['Model-1'].sketches['__profile__']
        mdb.models['Model-1'].parts.changeKey(fromName='Part-2', toName='BOX')

    def cut(self,part_name,z_min,mixPartName):
        a = mdb.models['Model-1'].rootAssembly
        a.DatumCsysByDefault(CARTESIAN)
        p = mdb.models['Model-1'].parts['BOX']
        a.Instance(name='BOX-1', part=p, dependent=ON)
        a = mdb.models['Model-1'].rootAssembly
        p = mdb.models['Model-1'].parts[part_name]
        instance_name=part_name+'-1'
        a.Instance(name=instance_name, part=p, dependent=ON)

        a = mdb.models['Model-1'].rootAssembly
        a.translate(instanceList=('BOX-1', ), vector=(0.0, 0.0, z_min))
        a = mdb.models['Model-1'].rootAssembly
        a.InstanceFromBooleanCut(name='Part-1', 
            instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances['BOX-1'], 
            cuttingInstances=(a.instances[instance_name], ), originalInstances=SUPPRESS)
        mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='Part-1-1', 
            toName='boxSubBcc')
        a = mdb.models['Model-1'].rootAssembly
        a.features[instance_name].resume()
        a1 = mdb.models['Model-1'].rootAssembly
        
        a1.InstanceFromBooleanMerge(name=mixPartName, instances=(
            a1.instances['boxSubBcc'], a1.instances[instance_name], ), 
            keepIntersections=ON, originalInstances=SUPPRESS, domain=GEOMETRY)

    def material(self,mixPartName):
        from material import createMaterialFromDataString
        createMaterialFromDataString('Model-1', 'Alsimg-thermo', '2023', 
            """{'description': '', 'density': {'temperatureDependency': OFF, 'table': ((2.7e-09,),), 'dependencies': 0, 'fieldName': '', 'distributionType': UNIFORM}, 'specificHeat': {'temperatureDependency': OFF, 'table': ((900000000.0,),), 'dependencies': 0, 'law': CONSTANTVOLUME}, 'materialIdentifier': '', 'conductivity': {'temperatureDependency': OFF, 'table': ((180.0,),), 'dependencies': 0, 'type': ISOTROPIC}, 'name': 'Alsimg-thermo'}""")
        from material import createMaterialFromDataString
        createMaterialFromDataString('Model-1', 'paraffin-solid', '2023', 
            """{'description': '', 'density': {'temperatureDependency': OFF, 'table': ((9e-10,),), 'dependencies': 0, 'fieldName': '', 'distributionType': UNIFORM}, 'specificHeat': {'temperatureDependency': OFF, 'table': ((2500000000.0,),), 'dependencies': 0, 'law': CONSTANTVOLUME}, 'materialIdentifier': '', 'latentHeat': {'table': ((200.0, 57.0, 63.0),)}, 'conductivity': {'temperatureDependency': OFF, 'table': ((0.3,),), 'dependencies': 0, 'type': ISOTROPIC}, 'name': 'paraffin-solid'}""")
        del mdb.models['Model-1'].materials['paraffin-solid'].latentHeat
        mdb.models['Model-1'].materials['paraffin-solid'].specificHeat.setValues(
            temperatureDependency=ON, table=((2500000000.0, 25.0), (2500000000.0, 
            54.0), (22500000000.0, 55.0), (22500000000.0, 65.0), (2500000000.0, 
            66.0), (2500000000.0, 90.0)))
        mdb.models['Model-1'].HomogeneousSolidSection(name='Al', 
            material='Alsimg-thermo', thickness=None)
        p1 = mdb.models['Model-1'].parts[mixPartName]
        mdb.models['Model-1'].HomogeneousSolidSection(name='Para', 
            material='paraffin-solid', thickness=None)
        p = mdb.models['Model-1'].parts[mixPartName]
        c = p.cells
        cells = c.getSequenceFromMask(mask=('[#2 ]', ), )
        region = p.Set(cells=cells, name='Al')
        p = mdb.models['Model-1'].parts[mixPartName]
        p.SectionAssignment(region=region, sectionName='Al', offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', 
            thicknessAssignment=FROM_SECTION)
        p = mdb.models['Model-1'].parts[mixPartName]
        c = p.cells
        cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        region = p.Set(cells=cells, name='Para')
        p = mdb.models['Model-1'].parts[mixPartName]
        p.SectionAssignment(region=region, sectionName='Para', offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', 
            thicknessAssignment=FROM_SECTION)

stp_path='C:/Users/Administrator/Desktop/20251022/calMeltingTime/137777.stp'
part_name='originStructure'
mixPartName='mix'
stpImport=stp2mixInstance(stp_path,part_name,mixPartName)
x_max,y_max,z_max,x_middle,y_middle,z_middle,x_min,y_min,z_min=stpImport.importStp(stp_path,part_name)
stpImport.buildBox(x_max,y_max,z_max,x_min,y_min,z_min)
stpImport.cut(part_name,z_min,mixPartName)
stpImport.material(mixPartName)


class generateSymmetryMesh:
    def __init__(self,part_name,output_name):
        self.part_name = part_name    
        self.job_name = os.path.splitext(os.path.basename(output_name))[0]
    def importstp(self,part_name):
        print(part_name)
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

    def buildbox(self,x_max,x_middle,y_max,y_middle,z_max,z_middle):

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

    def cut(self,part_name,z_middle):
        a = mdb.models['Model-1'].rootAssembly
        a.DatumCsysByDefault(CARTESIAN)
        p = mdb.models['Model-1'].parts['box']
        a.Instance(name='box-1', part=p, dependent=ON)
        p = mdb.models['Model-1'].parts[part_name]
        instance_name=part_name+'-1'
        a.Instance(name=instance_name, part=p, dependent=ON)
        a.translate(instanceList=('box-1', ), vector=(0.0, 0.0, z_middle))

        print('geometryValidity',a.instances[instance_name].part.geometryValidity) 

        a.InstanceFromBooleanCut(name='Part-boxSub8th', 
            instanceToBeCut=a.instances[instance_name], 
            cuttingInstances=(a.instances['box-1'], ), originalInstances=SUPPRESS)
        
        a.features[instance_name].resume()
        a.InstanceFromBooleanCut(name='Part-8th', 
            instanceToBeCut=a.instances[instance_name], 
            cuttingInstances=(a.instances['Part-boxSub8th-1'], ), 
            originalInstances=SUPPRESS)
        #p = mdb.models['Model-1'].parts['e1']
    
    def wangge(self):
        
        p = mdb.models['Model-1'].parts['Part-8th']
        c = p.cells
        # pickedRegions = c.getSequenceFromMask(mask=('[#1 ]', ), )
        pickedRegions = c
        p.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
        elemType3 = mesh.ElemType(elemCode=C3D10)
        p = mdb.models['Model-1'].parts['Part-8th']
        c = p.cells
        # cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        cells = c
        pickedRegions =(cells, )
        p.setElementType(regions=pickedRegions, elemTypes=( 
            elemType3,))
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
                  
    def mirror(self,x_middle,y_middle,z_middle):
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
        a.features['Part-8th-MirrorXY-1'].suppress()
        # a = mdb.models['Model-1'].rootAssembly
        a.InstanceFromBooleanMerge(name='mirror-finished', instances=(a.instances['Part-2nd-1'], 
            a.instances['Part-2nd-copy-MirrorXY-1'], ), mergeNodes=BOUNDARY_ONLY, 
            nodeMergingTolerance=1e-06, domain=MESH, originalInstances=SUPPRESS)

    def inpex(self):
        job_name=self.job_name

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



output_name = os.getenv("BPMESH_OUTPUT_MESH", r'C:\Users\Administrator\Desktop\20251022\calMeltingTime\137777-mixSymmetryMesh.inp') #output inp
symmetryMesh=generateSymmetryMesh(mixPartName,output_name)
x_max,y_max,z_max,x_middle,y_middle,z_middle=symmetryMesh.importstp(mixPartName)
x_center=x_middle
y_center=y_middle
z_center=z_middle

symmetryMesh.buildbox(x_max,x_center,y_max,y_center,z_max,z_center)
symmetryMesh.cut(mixPartName,z_middle)
symmetryMesh.wangge()
symmetryMesh.mirror(x_center,y_center,z_center)
symmetryMesh.inpex()
print(output_name+' finished!')