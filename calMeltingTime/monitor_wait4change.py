#!/usr/bin/env python
# coding=utf-8

"""
This script is used to monitor abaqus odb result during the simulation

EXIT CODE:
    - 1 : not ready and no error, keep waiting, and export csv also
    - 0 : Stop computing due to early stop criteria or errorneous state


Maintainer: Zhang Junxi, Zhang Ruixuan

-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
2025-09-23  Add check for non-physical issues
2025-09-10  init
"""

import sys
import os
import csv
import json
import time
import numpy as np

import odbAccess
from abaqusConstants import *


def get_RF_series():
    """
    Get reaction force series of corresponding boundary
    """
    print("")
    print("")
    print("[%s] >>> get_RF_series()" % time.ctime())

    xy2 = []

    for i_step, step in enumerate(odb.steps.values()):
        frames = step.frames
        for i_frame, frame in enumerate(frames):
            if i_step > 0 and i_frame == 0:
                continue  # skip the first frame of each step except the first step
            rfField = frame.fieldOutputs['RF'].getSubset(region=nodeSet_DISP, position=NODAL)
            assert len(rfField.bulkDataBlocks) == 1, "More than 1 block found in RF field output!"
            rf_sum = rfField.bulkDataBlocks[0].data[:, idx].sum()
            xy2.append((frame.frameValue + (0 if i_step == 0 else 1), abs(rf_sum)))
    
    if not xy2:
        print("No data found in odb!")
        sys.exit(1)
    
    x2, y2 = zip(*xy2)
    x2 = list(x2)
    y2 = list(y2)
    print("[%s] steps=%d, last time=%f" % (time.ctime(), len(xy2), x2[-1]))
    return x2, y2

def get_U_last():
    """
    2025-09-24  Bugfix: Step-2 always exists if set, but may have 0 frames
    """
    print("[%s] >>> get_U_last()" % time.ctime())

    step = odb.steps['Step-1']
    frameAcc = 0
    if 'Step-2' in odb.steps:
        if len(odb.steps['Step-2'].frames) > 0:
            step = odb.steps['Step-2']
            frameAcc = 1
    frame = step.frames[-1]
    u_field = frame.fieldOutputs['U'].getSubset(region=nodeSet_BOUNDARY, position=NODAL)
    assert len(u_field.bulkDataBlocks) == 1, "More than 1 block found in U field output!"
    U = u_field.bulkDataBlocks[0].data
    return U, frame.frameValue + frameAcc


def get_coords_boundary():
    print("[%s] >>> get_coords_boundary()" % time.ctime())
    coords = np.array([n.coordinates for n in nodeSet_BOUNDARY.nodes[0]])
    assert coords.shape[1] == 3, "Coordinates should be 3D!"
    return coords


def IsUValid(coords, U):
    """
    1. check if any node moves outside bounding box
    2. check if disp-boundary lower than others
    3. check if any self-collision happens
    """
    print("[%s] >>> IsUValid()" % time.ctime())

    newU = coords + U
    maxMD_DISP = newU[DISP_in_BOUNDARY, idx].max()
    maxMD_BOUNDARY = newU[:, idx].max()
    if maxMD_BOUNDARY > maxMD_DISP + 2e-1:
        print("[%s] Non-physical issue detected: max disp of boundary nodes is larger than disp nodes!" % (time.ctime()))
        print("[%s]   - maxMD_DISP = %.3f" % (time.ctime(), maxMD_DISP))
        print("[%s]   - maxMD_BOUNDARY = %.3f" % (time.ctime(), maxMD_BOUNDARY))
        return False 

    if bcmode > 0:
        TopBot_idxs = np.any(np.isclose(np.abs(newU), 5), axis=1)
        NoTopBot_idxs = np.setdiff1d(np.arange(newU.shape[0]), TopBot_idxs)
        if np.any(np.abs(newU[NoTopBot_idxs, :]) > 5 + 1e-2):
            print("[%s] Non-physical issue detected: some internal nodes moved outside the bounding box!" % (time.ctime()))
            print("[%s] Goes to: %.3f" % (time.ctime(), np.abs(newU[NoTopBot_idxs, :]).max()))
            return False

    if not os.path.exists("boundaryNodes.npy"):
        np.save("boundaryNodes.npy", nodeIds_BOUNDARY)
    np.save("newU.npy", newU)
    print("[%s] Gonna run: %s %s" % (time.ctime(), MSHPY, SCRIPT_CSC))
    istat = os.system("unset PYTHONPATH && %s %s" % (MSHPY, SCRIPT_CSC))
    if istat != 0:
        print("[%s] Non-physical issue detected: self-collision happens!" % (time.ctime()))
        return False
    
    return True


def export_csv(xs, ys):
    disp_data = []
    vflags = []
    for x in xs:
        if x < 1:
            disp_data.append(x * max_disp1)
        else:
            assert max_disp2 is not None, "max_disp2 should be set if time goes beyond 1!"
            disp_data.append(max_disp1 + (x - 1) * (max_disp2 - max_disp1))
        if x < abs(info["validStep"]) + 1e-5:
            vflags.append(1)
        else:
            vflags.append(0)
    disp_data_positive = [abs(x) for x in disp_data]

    combined_data = zip(disp_data_positive, ys, vflags)

    with open(outcsv, mode='w') as csv_file:
        writer = csv.writer(csv_file)
        if info:
            csv_file.write("% n_element=" + str(info["n_element"]) + ",n_node=" + str(info["n_node"]) + ",solve_time=" + str(time_cost) + ",cpus=" + str(info["ncpus"]) + ",validStep=" + str(abs(info["validStep"])) + "\n")
        writer.writerow(["-Displacement_U" + direction + "(mm)", "Reactionforce_RF" + direction + "(N)", "ValidFlag"])
        writer.writerows(combined_data)


direction = os.environ.get("PLASIM_DIRECTION", "x")
max_disp_env = os.environ.get("PLASIM_MAX_DISP", "-0.3,-1")
bcmode = int(os.environ.get("PLASIM_BCMODE", "0"))  # 0 means no constraint except for main axis, 1 means normal constraint at side faces; 2 means tagent constraint at all faces
MSHPY = os.environ.get("MSHPY", "/public/home/zncszjx/Software/idep/miniconda3/envs/fenicsx/bin/python3.12")
SCRIPT_CSC = os.environ.get("PLASIM_SCSC", "/public/home/nieqi01/zrx/abaqus/nonlinear/Lattice-Elasto-Plastic-Sim/Abaqus/check_self_collision.py")

print("")
print("====================================================")
print("[%s] Start monitor.py, direction=%s, max_disp_env=%s, bcmode=%d" % (time.ctime(), direction, max_disp_env, bcmode))

if "," not in max_disp_env:
    max_disp1 = float(max_disp_env)
    max_disp2 = None
else:
    max_disp1, max_disp2 = map(float, max_disp_env.split(","))

idx = 0 if direction == "x" else (1 if direction == "y" else 2)
outcsv = "displacement_RF%d.csv" % (idx+1)
odbfile = "Job-1.odb"

if not os.path.exists(odbfile):
    print("[%s] %s not found!" % (time.ctime(), odbfile))
    sys.exit(1)


iCount = 0
while True:
    if iCount > 10:
        print("[%s] Odb file is being written for too long, STOP computation" % (time.ctime()))
        sys.exit(0)
    iCount += 1
    try:
        odb = odbAccess.openOdb(path=odbfile, readOnly=True)
        break
    except odbAccess.OdbError:
        print("[%s] Odb file is being written, wait for 10s ..." % (time.ctime()))
        time.sleep(5)

nodeSet_DISP = odb.rootAssembly.nodeSets[direction.upper() + "TOP"]
nodeSet_BOUNDARY = odb.rootAssembly.nodeSets["ALLBOUNDARY"]
nodeIds_DISP = np.array([n.label for n in nodeSet_DISP.nodes[0]]) - 1
nodeIds_BOUNDARY = np.array([n.label for n in nodeSet_BOUNDARY.nodes[0]]) - 1

DISP_in_BOUNDARY = np.searchsorted(nodeIds_BOUNDARY, nodeIds_DISP)


info = None
if os.path.exists("info.json"):
    with open("info.json", "r") as f:
        info = json.load(f)
    time_cost = time.time() - info["start_time"]
assert info is not None


if info["validStep"] >= 0:
    coords = get_coords_boundary()
    U, frval = get_U_last()
    print("frval = ", frval)
    if not IsUValid(coords, U):
        print("[%s] Invalid displacement field detected, mark it" % (time.ctime()))
        info["validStep"] = -info["validStep"]
    else:
        info["validStep"] = frval
    with open("info.json", "w") as f:
        json.dump(info, f)

xs, ys = get_RF_series()
export_csv(xs, ys)

maxTime = 2 if max_disp2 is not None else 1


if abs(xs[-1] - maxTime) < 1e-5:
    print("[%s] Max time reached, STOP computation" % (time.ctime()))
    sys.exit()
elif len(ys) > 10 and ys[-1] < ys[-2] and ys[-2] < ys[-3] and ys[-3] < ys[-4]:
    print("[%s] 4 decreasing points detected, STOP computation" % (time.ctime()))
    sys.exit()
else:
    print("[%s] Still increasing (t[-1] = %f)" % (time.ctime(), xs[-1]))
    print("ys[-6:] = ", ys[-6:])
    print("[%s] continue ..." % (time.ctime()))
    sys.exit(1)
