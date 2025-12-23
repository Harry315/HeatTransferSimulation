from abaqus import *
from abaqusConstants import *
import xyPlot
from odbAccess import openOdb
import numpy as np


#Use Macro to get the value of the HFL2 node by node! 

def getHFL_gui_pc(odb_path,set_name):
    HFLname = (
    'HFL1' if 'X' in set_name else
    'HFL2' if 'Y' in set_name else
    'HFL3' if 'Z' in set_name else
    None
    )

    o1 = session.openOdb(
        name=odb_path)
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    odb = session.odbs[odb_path]
    session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('HFL', 
        INTEGRATION_POINT, ((COMPONENT, HFLname), )), ), nodeSets=(set_name, ))


    xy1 = session.xyDataObjects
    node_names=xy1.keys()  
    node_num=len(node_names)
    HFL_sum=0
    for node_name in node_names:
        HFL_node=xy1[node_name].data[-1][-1]
        HFL_sum+=HFL_node
        del xy1[node_name]

    print('average HFL gui',HFL_sum/node_num)
    print('point number:',node_num)
    #print(HFL_sum/node_num*2)
    

def getHFL_gui(odb_path,set_name):
    HFLname = (
    'HFL1' if 'X' in set_name else
    'HFL2' if 'Y' in set_name else
    'HFL3' if 'Z' in set_name else
    None
    )

    odb = session.odbs[odb_path]

    xy2=xyPlot.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('HFL', 
        INTEGRATION_POINT, ((COMPONENT, HFLname), )), ), nodeSets=(set_name, ))
    node_num=len(xy2)
    HFL_sum=0
    for i in range(node_num):
        HFL_node=xy2[i].data[1][1]
        HFL_sum+=HFL_node


    print('average HFL gui',HFL_sum/node_num)
    print('point number:',node_num)

    
def getHFL_nogui(odb_path,set_name):
    odb = openOdb(path=odb_path, readOnly=True)

    last_step = odb.steps[odb.steps.keys()[-1]]
    last_frame = last_step.frames[-1]

    node_set = odb.rootAssembly.nodeSets[set_name]
    print('node_set',len(node_set.nodes[0]))


    node_labels = [n.label for n in node_set.nodes[0]]  


    hfl_nodal = last_frame.fieldOutputs['HFL'].getSubset(position=ELEMENT_NODAL)
    hfl_top   = hfl_nodal.getSubset(region=node_set)

    last_val = {}
    for v in hfl_top.values:
        last_val[v.nodeLabel] = v.data[1]  

    hfl_array = np.array([last_val[lbl] for lbl in node_labels if lbl in last_val])
    print('average HFL nogui:',hfl_array.mean())
    print('point number:',len(hfl_array))

    odb.close()

def getHFL_nogui_weight(odb_path):
    odb = openOdb(path=odb_path, readOnly=True)

    last_step = odb.steps[odb.steps.keys()[-1]]
    frame = last_step.frames[-1]

    node_set = odb.rootAssembly.nodeSets['XTOP']
    target_lbl = {n.label for n in node_set.nodes[0]}   

    elemMap = {}
    for inst in odb.rootAssembly.instances.values():
        for e in inst.elements:
            elemMap[e.label] = (inst, e)

    hfl_nodal = frame.fieldOutputs['HFL'].getSubset(position=ELEMENT_NODAL)
    hfl_top   = hfl_nodal.getSubset(region=node_set)

    sumVal = {}   
    sumVol = {}   

    for v in hfl_top.values:
        if v.elementLabel is None:
            continue    
        nid = v.nodeLabel
        if nid not in target_lbl:
            continue
        inst, elem = elemMap[v.elementLabel]
        vol = elem.volume
        w   = vol / len(elem.connectivity)   
        sumVal[nid] = sumVal.get(nid, 0.0) + v.data[0] * w
        sumVol[nid] = sumVol.get(nid, 0.0) + w


    weighted = np.array([sumVal[lbl]/sumVol[lbl] for lbl in sumVal])
    print('average HFL nogui_weight:',weighted.mean())  
    print('point number:',len(weighted))

    odb.close()

def HFL_weight_nogui(odb_path,set_name):
    #First, obtain the nodes of each surface element and compare them with the existing node list. 
    # Since some elements reuse the same node, this increases the weight.
    
    odb = openOdb(path=odb_path, readOnly=True)       
    node_set = odb.rootAssembly.nodeSets[set_name]
    print('node_set',len(node_set.nodes[0]))
    node_labels_weight = [n.label for n in node_set.nodes[0]]  #all nodes
    # print(node_labels_weight)
    nodeWeight = {label: 0 for label in node_labels_weight}

    surfFaces = odb.rootAssembly.instances['MIRROR-FINISHED-1'].surfaces['SURF-1']
    lelements=len(surfFaces.elements)
    for i in range(lelements):
        # print(surfFaces.elements[i].connectivity)
        for nd in surfFaces.elements[i].connectivity:
            if nd in nodeWeight:          
                nodeWeight[nd] += 1
    # print(nodeWeight)     

    last_step  = odb.steps[odb.steps.keys()[-1]]
    last_frame = last_step.frames[-1]

    node_set  = odb.rootAssembly.nodeSets[set_name]
    node_labels = [n.label for n in node_set.nodes[0]]   

    hfl_nodal = last_frame.fieldOutputs['HFL'].getSubset(position=ELEMENT_NODAL)
    hfl_top   = hfl_nodal.getSubset(region=node_set)


    hfl_val = {} 
    for v in hfl_top.values:
        hfl_val[v.nodeLabel] = v.data[1]

    weight_sum = 0.0
    hfl_weighted_sum = 0.0
    for lbl in node_labels:
        if lbl in hfl_val and lbl in nodeWeight:
            w = nodeWeight[lbl]
            weight_sum += w
            hfl_weighted_sum += w * hfl_val[lbl]

    if weight_sum == 0:
        avg_hfl = np.nan
    else:
        avg_hfl = hfl_weighted_sum / weight_sum

    print('nogui_HFL_surface:', avg_hfl)
    print('weight_sum:', weight_sum)
    odb.close()

def HFL_weight_gui(odb_path,set_name):

    odb = openOdb(path=odb_path, readOnly=True)       
    node_set = odb.rootAssembly.nodeSets[set_name]
    print('node_set',len(node_set.nodes[0]))
    node_labels_weight = [n.label for n in node_set.nodes[0]]  #all nodes
    # print(node_labels_weight)
    nodeWeight = {label: 0 for label in node_labels_weight}

    surfFaces = odb.rootAssembly.instances['MIRROR-FINISHED-1'].surfaces['SURF-1']
    lelements=len(surfFaces.elements)
    for i in range(lelements):
        # print(surfFaces.elements[i].connectivity)
        for nd in surfFaces.elements[i].connectivity:
            if nd in nodeWeight:          
                nodeWeight[nd] += 1
    # print(nodeWeight)     

    # weighted average HFL2 based on nodeWeight
    o1 = session.openOdb(name=odb_path)
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    odb = session.odbs[odb_path]

    session.xyDataListFromField(
        odb=odb,
        outputPosition=NODAL,
        variable=(('HFL', INTEGRATION_POINT, ((COMPONENT, 'HFL2'),)),),
        nodeSets=(set_name,)
    )

    xy1 = session.xyDataObjects
    node_names = xy1.keys()

    weight_sum = 0.0
    hfl_weighted = 0.0

    for curve_name in node_names:
        nodeLabel = int(curve_name.split()[-1])
        if nodeLabel not in nodeWeight:
            continue

        w = nodeWeight[nodeLabel]
        hfl_val = xy1[curve_name].data[-1][-1]

        weight_sum += w
        hfl_weighted += w * hfl_val
        del xy1[curve_name]

    avg_hfl = hfl_weighted / weight_sum if weight_sum else float('nan')
    print('Weighted average HFL:', avg_hfl)
    print('Total weight:', weight_sum)




#odb_path = r'C:/Users/Administrator/Desktop/20251022/abaqus_heat/box-stationary-25-30.odb'
odb_path = r'C:/Users/Administrator/Desktop/20251022/abaqus_heat/stationary-25-30-symmetryMesh.odb'
odb_path = r'C:\Users\Administrator\Desktop\20251022\paraffin_mix\bcc-mix-symmetryMesh\mix-symmetryMesh.odb'
odb_path = r'C:\Users\Administrator\Desktop\20251022\paraffin_mix\bcc-mix-symmetryMesh\mix-thermo-DC3D4.odb'
odb_path = r'E:\workspace\20251022\calHFL\137777-symmetryMesh.odb'
#odb_path = r'C:\Users\Administrator\Desktop\20251022\paraffin_mix\bcc-mix-symmetryMesh\mixChangeProperty.odb'
#odb_path = r'C:\Users\Administrator\Desktop\20251022\paraffin_mix\bcc-mix-symmetryMesh\mixTemp25-50.odb'
#odb_path = r'C:/Users/Administrator/Desktop/20251022/abaqus_heat/staionary-25-30-array222.odb'
print(odb_path)

set_name='XTOP'
# set_name='XBOT'
print(set_name)
# print('HFL2')
getHFL_gui(odb_path,set_name)
# getHFL_nogui(odb_path,set_name)
# HFL_weight_nogui(odb_path,set_name)
# HFL_weight_gui(odb_path,set_name)
# getHFL_nogui_weight(odb_path)