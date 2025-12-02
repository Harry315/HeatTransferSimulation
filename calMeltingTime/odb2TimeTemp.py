"""
Extract average NT11 vs time for node-set 'YBOT' across Step-1 and Step-2
(step-2 continues from step-1)
"""

from abaqus import *
from abaqusConstants import *
import numpy as np

def processNodeSets():
    odb_path = 'mix-thermo-oneSide.odb'
    odb_path = 'mix-thermo-oneSideHFL.odb'
    odb_path = 'mix-thermo-oneSide-changeSpecificHeat.odb'
    #odb_path = 'mix-thermo-oneSideHFL-bigLatent.odb'
    # set_name = 'YBOT'
    set_name = 'ALL'
    print(odb_path)
    odb = session.openOdb(name=odb_path, readOnly=True)
    ra  = odb.rootAssembly

    avg_list  = []
    time_list = []

    for step_name in ['Step-1', 'Step-2']:
    #for step_name in ['Step-1']:
        step = odb.steps[step_name]
        for frame in step.frames:
            time_list.append(frame.frameValue)
            nt_field = frame.fieldOutputs['NT11']
            nt_set   = nt_field.getSubset(region=ra.nodeSets[set_name], position=NODAL)
            avg_list.append(np.mean([v.data for v in nt_set.values]))
    return avg_list,time_list
    #odb.close()



#====================================================

def processInstNodeSets(odb_path,set_name):

    inst_name, set_name = set_name.rsplit('.', 1)
    print(odb_path)
    odb = session.openOdb(name=odb_path, readOnly=True)
    ra  = odb.rootAssembly

    avg_list  = []
    time_list = []

    for step_name in ['Step-1', 'Step-2']:
    #for step_name in ['Step-1']:
        step = odb.steps[step_name]
        for frame in step.frames:
            time_list.append(frame.frameValue)
            nt_field = frame.fieldOutputs['NT11']
            nt_set   = nt_field.getSubset(region=ra.instances[inst_name].nodeSets[set_name], position=NODAL)
            avg_list.append(np.mean([v.data for v in nt_set.values]))
    return avg_list,time_list
    #odb.close()


#=================================================
odb_path = 'mix-thermo-oneSide-changeSpecificHeat.odb'

set_name = 'MIRROR-FINISHED-1.AL'
# set_name = 'MIRROR-FINISHED-1.PARA'
avg_list,time_list=processInstNodeSets(odb_path,set_name)


time_arr = np.array(time_list)
avg_arr  = np.array(avg_list)

print('Time (s):', time_arr)
print('Avg NT11 :', avg_arr)
print('\n----- For Excel -----\n')
for t, v in zip(time_arr, avg_arr):
    print(str(t)+' '+str(v))