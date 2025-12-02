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

def importInp(inp_path):

    filename_with_ext = os.path.basename(inp_path)

    filename_without_ext = os.path.splitext(filename_with_ext)[0]
    mdb.ModelFromInputFile(name=filename_without_ext, 
        inputFileName=inp_path)
    return filename_without_ext

def changeMesh(model_name,part_name):
    p = mdb.models[model_name].parts[part_name]
    elemType1 = mesh.ElemType(elemCode=DC3D10, elemLibrary=STANDARD)
    z1 = p.elements
    pickedRegions =(z1, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, ))

def set_boundary(model_name,instance_name,meshRange=1e-4):
    a = mdb.models[model_name].rootAssembly
    s1 = a.instances[instance_name].nodes
    x_vals = [n.coordinates[0] for n in s1]
    y_vals = [n.coordinates[1] for n in s1]
    z_vals = [n.coordinates[2] for n in s1]

    # [x_min, x_max, y_min, y_max, z_min, z_max]
    bounds = [min(x_vals), max(x_vals),
            min(y_vals), max(y_vals),
            min(z_vals), max(z_vals)]

    boundarySets=['xtop','ytop','ztop','xbot','ybot','zbot', 'allBoundary'] # ,'xhalf','yhalf','zhalf']
    for boundaryName in boundarySets:
        # a = mdb.models[MODELNAME].rootAssembly
        # s1 = a.instances[INSTNAME].faces
        coor = {}.fromkeys(['xMin', 'yMin', 'zMin', 'xMax', 'yMax', 'zMax'])
        coor['xMin'], coor['xMax'], coor['yMin'], coor['yMax'], coor['zMin'], coor['zMax']  = bounds
        
        boundaryDirection=boundaryName[0]
        if 'top' in boundaryName:
            coor['%sMin' % (boundaryDirection)]=coor['%sMax' % (boundaryDirection)]-meshRange
            coor['%sMax' % (boundaryDirection)]=coor['%sMax' % (boundaryDirection)]+meshRange
        elif 'bot' in boundaryName:
            coor['%sMax' % (boundaryDirection)]=coor['%sMin' % (boundaryDirection)]+meshRange
            coor['%sMin' % (boundaryDirection)]=coor['%sMin' % (boundaryDirection)]-meshRange

        faces_nodes = a.instances[instance_name].sets['AL'].nodes.getByBoundingBox (**coor)
        if len(faces_nodes)==0:
            print('No point has been selected!')
        print("Get %s: %d faces"%(boundaryName, len(faces_nodes)))
        a.Set(nodes=faces_nodes, name=boundaryName)

def tempBoundaryCondition(model_name,direction,instance_name):
    mdb.models[model_name].setValues(absoluteZero=-273.15)
    a = mdb.models[model_name].rootAssembly
    a.regenerate()
    mdb.models[model_name].HeatTransferStep(name='Step-1', 
        previous='Initial', timePeriod=0.001, initialInc=0.001, minInc=1e-08, 
        maxInc=0.001, deltmx=20.0, amplitude=RAMP)
    mdb.models[model_name].HeatTransferStep(name='Step-2', 
        previous='Step-1', timePeriod=400.0, initialInc=0.01, minInc=1e-08, 
        maxInc=400.0, deltmx=20.0)
    
    region = a.sets[direction+'top']
    mdb.models[model_name].TemperatureBC(name='BC-1', 
        createStepName='Step-1', region=region, fixed=OFF, 
        distributionType=UNIFORM, fieldName='', magnitude=65.0, 
        amplitude=UNSET)
    n1 = a.instances[instance_name].nodes
    region = a.Set(nodes=n1, name='all')
    mdb.models[model_name].Temperature(name='Predefined Field-1', 
        createStepName='Initial', region=region, distributionType=UNIFORM, 
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(25.0, 
        ))
    
def submit(model_name,ncpus=1):
    print("%s Solving..." % time.ctime())
    mdb.Job(
        name=model_name, model=model_name, description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=ncpus,
        numDomains=ncpus, numGPUs=0)
    mdb.jobs[model_name].submit(consistencyChecking=OFF)
    mdb.jobs[model_name].waitForCompletion()

infile=r'C:\Users\Administrator\Desktop\20251022\calMeltingTime\137777-mixSymmetryMesh.inp'
model_name=importInp(infile)
instance_name='MIRROR-FINISHED-1'
part_name='MIRROR-FINISHED'
direction="y"
ncpus=1

assert ncpus > 0
assert os.path.exists(infile)

changeMesh(model_name,part_name)
set_boundary(model_name,instance_name)
tempBoundaryCondition(model_name,direction,instance_name)
submit(model_name,ncpus)