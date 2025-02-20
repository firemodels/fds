# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 18:21:38 2024

@author: jhodges
"""
import numpy as np
import os
from collections import defaultdict

def getChid(testnum, convmodel, dx):
    if convmodel == 'jet':
        suffix1 = '_jet'
        suffix3 = '_%03dmm'%(dx*1e3)
    else:
        suffix1 = '_default'
        suffix3 = '_%03dmm'%(dx*1e3)

    chid = 'lattimer_test%0d%s%s'%(testnum, suffix1, suffix3)
    return chid

def splitMeshOnce(meshes, mesh, meshSplitAxes=[True, True, True]):
    IJK = np.round(mesh['IJK'])
    XB = mesh['XB']
    dxs = [(XB[1]-XB[0])/float(IJK[0]), (XB[3]-XB[2])/float(IJK[1]), (XB[5]-XB[4])/float(IJK[2])]
    ind = np.argmax(IJK)
    IJK_temp = list(IJK)
    while meshSplitAxes[ind] is False:
        IJK_temp = list(IJK_temp)
        IJK_temp[ind] = 0
        ind = np.argmax(IJK_temp)
        if np.sum(IJK_temp) == 0:
            print("Failed to split mesh.")
            break
    
    IJK2 = list(IJK)
    XB2 = list(XB)
    IJK2[ind] = int(IJK[ind]/2)
    if IJK2[ind] % 2 > 0: IJK2[ind] = IJK2[ind]-1
    XB2[int(2*ind+1)] = XB2[int(2*ind)] + dxs[ind]*float(IJK2[ind])
    
    IJK3 = list(IJK)
    XB3 = list(XB)
    IJK3[ind] = IJK[ind] - IJK2[ind]
    XB3[int(2*ind)] = XB2[int(2*ind+1)]
    
    mesh2 = defaultdict(bool)
    mesh2['ID'] = "%s-00"%(mesh["ID"])
    mesh2['XB'] = np.array(XB2)
    mesh2['IJK'] = np.array(IJK2)
    
    mesh3 = defaultdict(bool)
    mesh3['ID'] = "%s-01"%(mesh["ID"])
    mesh3['XB'] = np.array(XB3)
    mesh3['IJK'] = np.array(IJK3)
    
    meshes.pop(mesh['ID'], False)
    meshes[mesh2['ID']] = mesh2
    meshes[mesh3['ID']] = mesh3
    
    return meshes

def calculateMeshCells(meshes):
    mesh_ids = []
    numCells = []
    meshKeys = list(meshes.keys())
    try:
        meshKeys.remove('unknownCounter')
    except:
        pass
    for key in meshKeys:
        IJK = meshes[key]['IJK']
        numCells.append(IJK[0]*IJK[1]*IJK[2])
        mesh_ids.append(key)
    return mesh_ids, numCells

def buildMesh(XB, IJK, numberOfProcesses, splitMultiplier=1.2):
    meshes = defaultdict(bool)
    meshes['start'] = defaultdict(bool, {'ID': 'start', 'XB': XB, 'IJK': IJK})
    if numberOfProcesses == 1: return meshes
    cellsPerProcess = IJK[0]*IJK[1]*IJK[2]/numberOfProcesses
    splitConverged = False
    while not splitConverged:
        splitConverged = True
        mesh_ids, numCells = calculateMeshCells(meshes)
        for mesh, numCell in zip(mesh_ids, numCells):
            if numCell > cellsPerProcess*splitMultiplier:
                meshes = splitMeshOnce(meshes, meshes[mesh])
                splitConverged = False
        if len(meshes) / 10 > numberOfProcesses:
            print("Warning: Number of meshes 10x greater than number of requested processes (%0.0f, %0.0f)"%(len(meshes), numberOfProcesses))
            print("AssumingConvergence")
            splitConverged = True
    
    mesh_ids, numCells = calculateMeshCells(meshes)
    return meshes

def writeMeshes(meshes):
    mesh_ids, numCells = calculateMeshCells(meshes)
    txt = ''
    for i in range(0, len(mesh_ids)):
        IJK = meshes[mesh_ids[i]]['IJK']
        XB = meshes[mesh_ids[i]]['XB']
        IJK_s = '%d,%d,%d'%(IJK[0],IJK[1],IJK[2])
        XB_s = '%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,'%(XB[0],XB[1],XB[2],XB[3],XB[4],XB[5])
        ID = meshes[mesh_ids[i]]['ID']
        txt = txt + "&MESH ID='%s', IJK=%s, XB=%s /\n"%(ID, IJK_s, XB_s)
    return txt

def buildInputFile(testnum, HRR, H, convmodel, dx):
    import numpy as np
    chid = getChid(testnum, convmodel, dx)
    if convmodel == 'jet':
        suffix2 = "HEAT_TRANSFER_MODEL='IMPINGING JET', CONVECTION_LENGTH_SCALE=0.77, "
    else:
        suffix2 = ""
        
    txt = "&HEAD CHID='%s', /\n"%(chid)
    txt = txt + "&TIME T_END=200. /\n"
    txt = txt + "&DUMP DT_BNDF=1, DT_DEVC=0.1, WRITE_XYZ=T, RESULTS_DIR='results_%s', SUPPRESS_DIAGNOSTICS=T /\n"%(chid)
    txt = txt + "&RADI NUMBER_RADIATION_ANGLES=400, TIME_STEP_INCREMENT=12,  /\n"
    txt = txt + "&COMB SUPPRESSION=F /\n"
    txt = txt + "&MISC SIMULATION_MODE='LES' /\n"
    
    x = 3.2
    y = 1.8
    dz = dx
    H_MIN = 0.0
    H_MAX = 2.4
    
    K = (H_MAX-H_MIN)/dz
    I = np.round(x/dx)
    J = np.round(y/dx)
    
    if dx == 0.05: 
        numProcesses = 7
        coarse_dx = dx*1
        coarse_dz = dz*1
    if dx == 0.025: 
        numProcesses = 16
        coarse_dx = dx*2
        coarse_dz = dz*2
    if dx == 0.0125: 
        numProcesses = 32
        coarse_dx = dx*4
        coarse_dz = dz*4
    
    #H_MAX = np.round((K*dz+H_MIN+0.05)/coarse_dz)*coarse_dz
    XB = np.array([-x, 0.15, -y/2, y/2, H_MIN, H_MAX])
    I, J, K = np.round((XB[1]-XB[0])/dx), np.round((XB[3]-XB[2])/dx), np.round((XB[5]-XB[4])/dz)
    K = np.ceil(K/16)*16
    H_MAX = H_MIN+K*dz
    XB[5] = H_MAX
    IJK = np.array([I, J, K])
    
    Nx = np.ceil(I/7)
    xmn = -x
    for ii in range(0, 7):
        xmx = xmn + Nx*dx
        if xmx > x/2:
            xmx = x/2
            Nx = (xmx-xmn)/dx
        XB = np.array([xmn, xmx, -y/2+dx/2, y/2+dx/2, H_MIN, H_MAX])
        IJK = np.array([Nx, J, K])
        meshes = buildMesh(XB, IJK, numProcesses/7, splitMultiplier=1.2)
        txt = txt + writeMeshes(meshes)
        xmn = xmx

    txt = txt + "&MATL ID= 'CONCRETE', FYI = 'Drysdale, Intro to Fire Dynamics', SPECIFIC_HEAT = 0.90, DENSITY = 2000., CONDUCTIVITY  = 1.0, /\n"
    txt = txt + "&MATL ID='FIBERFRAX DURABLANKET S', SPECIFIC_HEAT=1.14, DENSITY=96, CONDUCTIVITY=0.14 /\n"
    
    txt = txt + "&SURF ID = 'FLOOR', MATL_ID = 'CONCRETE', THICKNESS  = 0.1, COLOR = 'GRAY' /\n"
    txt = txt + "&REAC PRIORITY=1, FUEL='PROPANE PILOT', SOOT_YIELD=0.024, CO_YIELD=0.005, RADIATIVE_FRACTION=0.29, AUTO_IGNITION_TEMPERATURE=0 / \n"
    txt = txt + "&REAC PRIORITY=1, FUEL='PROPANE FUEL', SOOT_YIELD=0.024, CO_YIELD=0.005, RADIATIVE_FRACTION=0.29, AUTO_IGNITION_TEMPERATURE=450 / \n"
    txt = txt + "&SPEC ID='PROPANE', /\n"
    txt = txt + "&SPEC ID='PROPANE PILOT', SPEC_ID='PROPANE', MASS_FRACTION=1/\n"
    txt = txt + "&SPEC ID='PROPANE FUEL', SPEC_ID='PROPANE', MASS_FRACTION=1/\n"
    
    txt = txt + "&SURF ID='CEILING', EMISSIVITY=0.8, MATL_ID='FIBERFRAX DURABLANKET S', THICKNESS=0.0254, COLOR = 'GRAY 80', %s /\n"%(suffix2)
    txt = txt + "&SURF ID='WALL', EMISSIVITY=0.8, MATL_ID='FIBERFRAX DURABLANKET S', THICKNESS=0.0254, COLOR='GRAY 40' /\n"
    txt = txt + "&SURF ID='fire', HRRPUA=%0.4f, TMP_FRONT=400, COLOR='RED', SPEC_ID='PROPANE PILOT','PROPANE FUEL', MASS_FRACTION=0.01,0.99, MASS_FLUX_VAR=0.1/\n"%(HRR/(0.495))
    
    txt = txt + "&OBST ID='CEILING', XB=-2.4,%0.4f,%0.4f,%0.4f,2.1,%0.4f, SURF_IDS='WALL','WALL','CEILING', /\n"%(dx,-0.6-dx,0.6+dx,2.1+dz)
    txt = txt + "&OBST ID='BACK-WALL', XB=0.0,%0.4f,%0.4f,%0.4f,0.9,2.1, SURF_ID='WALL', /\n"%(0.0+dx, -0.6-dx, 0.6+dx)
    txt = txt + "&OBST ID='SIDEWALL-1', XB=-2.4,0.0,-0.6,%0.4f,1.5,2.1, SURF_ID='WALL', /\n"%(-0.6-dx)
    txt = txt + "&OBST ID='SIDEWALL-2', XB=-2.4,0.0,0.6,%0.4f,1.5,2.1, SURF_ID='WALL', /\n"%(0.6+dx)
    txt = txt + "&OBST ID='SIDEWALL-3', XB=-1.2,0.0,-0.6,%0.4f,0.9,1.5, SURF_ID='WALL', /\n"%(-0.6-dx)
    txt = txt + "&OBST ID='SIDEWALL-2', XB=-1.2,0.0,0.6,%0.4f,0.9,1.5, SURF_ID='WALL', /\n"%(0.6+dx)
    
    txt = txt + "&OBST ID='FIRE', XB=-0.45,0.0,-0.55,0.55,%0.4f,%0.4f, SURF_IDS='fire','WALL','WALL' /\n"%(2.1-H-0.15,2.1-H)
    txt = txt + "&OBST ID='BLOCK', XB=-0.45,0.0,-0.55,0.55,0.0,%0.4f, SURF_ID='FLOOR' /\n"%(2.1-H-0.15)
    
    '''
    txt = txt + "&DEVC ID='t', QUANTITY='TIME', XYZ=0,0,0 /\n"
    txt = txt + "&CTRL ID='2*pi*t/10', FUNCTION_TYPE='MULTIPLY', INPUT_ID='CONSTANT','t', CONSTANT=0.62831853 /\n"
    txt = txt + "&CTRL ID='sin(2*pi*t/10)', FUNCTION_TYPE='SIN', INPUT_ID='2*pi*t/10' /\n"
    txt = txt + "&RAMP ID='noise', T=-1, F=0.8, CTRL_ID='sin(2*pi*t/10)' /\n"
    txt = txt + "&RAMP ID='noise', T= 1, F= 1 /\n"
    p_ramp = '' if wind > 0 else ", PRESSURE_RAMP='noise' "
    txt = txt + "&VENT DB='YMIN', SURF_ID='OPEN', DYNAMIC_PRESSURE= %0.4f%s /\n"%(abs(wind), p_ramp)
    txt = txt + "&VENT DB='YMAX', SURF_ID='OPEN', DYNAMIC_PRESSURE=-%0.4f%s /\n"%(abs(wind), p_ramp)
    '''
    
    txt = txt + "&VENT DB='XMIN', SURF_ID='OPEN', /\n"
    txt = txt + "&VENT DB='XMAX', SURF_ID='OPEN', /\n"
    txt = txt + "&VENT DB='YMIN', SURF_ID='OPEN', /\n"
    txt = txt + "&VENT DB='YMAX', SURF_ID='OPEN', /\n"
    txt = txt + "&VENT DB='ZMIN', SURF_ID='FLOOR', /\n"
    txt = txt + "&VENT DB='ZMAX', SURF_ID='OPEN', /\n"
    
    for i, xx in enumerate([-0.3, -0.9, -1.5, -2.1]):
        txt = txt + "&DEVC ID='GHF%d', XYZ=%0.4f,0.0,2.1, IOR=-3, QUANTITY='GAUGE HEAT FLUX', /\n"%(i, xx)
        txt = txt + "&DEVC ID='IHF%d', XYZ=%0.4f,0.0,2.1, IOR=-3, QUANTITY='INCIDENT HEAT FLUX', /\n"%(i, xx)
        txt = txt + "&DEVC ID='CHF%d', XYZ=%0.4f,0.0,2.1, IOR=-3, QUANTITY='CONVECTIVE HEAT FLUX', /\n"%(i, xx)
        txt = txt + "&DEVC ID='GCHF%d', XYZ=%0.4f,0.0,2.1, IOR=-3, QUANTITY='CONVECTIVE HEAT FLUX GAUGE', /\n"%(i, xx)
        txt = txt + "&DEVC ID='HTC%d', XYZ=%0.4f,0.0,2.1, IOR=-3, QUANTITY='HEAT TRANSFER COEFFICIENT', /\n"%(i, xx)
        txt = txt + "&DEVC ID='TGAS%d', XYZ=%0.4f,0.0,2.1, IOR=-3, QUANTITY='GAS TEMPERATURE', /\n"%(i, xx)
        txt = txt + "&DEVC ID='TWALL%d', XYZ=%0.4f,0.0,2.1, IOR=-3, QUANTITY='WALL TEMPERATURE', /\n"%(i, xx)
        txt = txt + "&DEVC ID='VELOCITY%d', XYZ=%0.4f,0.0,2.05, QUANTITY='VELOCITY', /\n"%(i, xx)
        txt = txt + "&DEVC ID='VELOCITY%d', XYZ=%0.4f,0.0,2.075, QUANTITY='TEMPERATURE', /\n"%(i, xx)
        
    
    specs = ['PROPANE PILOT', 'PROPANE FUEL']
    spec_counter = 0
    for qty in ['HEAT TRANSFER COEFFICIENT', 'GAS TEMPERATURE', 'ADIABATIC SURFACE TEMPERATURE', 'GAUGE HEAT FLUX', 'INCIDENT HEAT FLUX', 'CONVECTIVE HEAT FLUX', 'WALL TEMPERATURE', 'CONVECTIVE HEAT TRANSFER REGIME']:
        txt = txt + "&BNDF QUANTITY='%s', "%(qty)
        if qty == 'MASS FLUX':
            txt = txt + "SPEC_ID='%s', "%(specs[spec_counter])
            spec_counter += 1
        txt = txt + "CELL_CENTERED=T /\n"
    for qty in ['TEMPERATURE', 'HRRPUV']:
        txt = txt + "&SLCF PBX=0, QUANTITY='%s', "%(qty)
        txt = txt + "CELL_CENTERED=T /\n"
    txt = txt + "&DEVC ID='RTE_C', XYZ=0.0,0.0,0.5, QUANTITY='RTE SOURCE CORRECTION FACTOR' /\n"
    txt = txt + '&TAIL /\n'
    
    return txt, numProcesses

def generateSlurm(namespace, casenum, numberOfProcesses, chids, queue, runtime):
    runName = '%s%d'%(namespace, casenum)
    totalProcesses = np.sum(numberOfProcesses)
    txt = '#!/bin/bash\n'
    txt = txt + '#----------------------------------------------------\n'
    txt = txt + '# Sample Slurm job script\n'
    txt = txt + '#   for TACC Frontera CLX nodes\n'
    txt = txt + '#\n'
    txt = txt + '#   *** Serial Job in Small Queue***\n'
    txt = txt + '#\n'
    txt = txt + '# Last revised: 22 June 2021\n'
    txt = txt + '#\n'
    txt = txt + '# Notes:\n'
    txt = txt + '#\n'
    txt = txt + '#  -- Copy/edit this script as desired.  Launch by executing\n'
    txt = txt + '#     "sbatch clx.serial.slurm" on a Frontera login node.\n'
    txt = txt + '#\n'
    txt = txt + '#  -- Serial codes run on a single node (upper case N = 1).\n'
    txt = txt + '#       A serial code ignores the value of lower case n,\n'
    txt = txt + '#       but slurm needs a plausible value to schedule the job.\n'
    txt = txt + '#\n'
    txt = txt + '#  -- Use TACCs launcher utility to run multiple serial\n'
    txt = txt + '#       executables at the same time, execute "module load launcher"\n'
    txt = txt + '#       followed by "module help launcher".\n'
    txt = txt + '#----------------------------------------------------\n'
    txt = txt + '\n'
    txt = txt + '#SBATCH -J %s    # Job name\n'%(runName)
    txt = txt + '#SBATCH -o %s.o%%j       # Name of stdout output file\n'%(runName)
    txt = txt + '#SBATCH -e %s.e%%j       # Name of stderr error file\n'%(runName)
    txt = txt + '#SBATCH -p %s         # Queue (partition) name\n'%(queue)
    txt = txt + '#SBATCH -N 2              # Total # of nodes (must be 1 for serial)\n'
    txt = txt + '#SBATCH -n %d             # Total # of mpi tasks (should be 1 for serial)\n'%(totalProcesses)
    txt = txt + '#SBATCH -t %s        # Run time (hh:mm:ss)\n'%(runtime)
    txt = txt + '#SBATCH --mail-type=all    # Send email at begin and end of job\n'
    txt = txt + '#SBATCH -A PHY23048      # Project/Allocation name (reqd if you have more than 1)\n'
    txt = txt + '#SBATCH --mail-user=jhodges@jensenhughes.com\n'
    txt = txt + '\n'
    txt = txt + '# Any other commands must follow all #SBATCH directives...\n'
    txt = txt + 'module load intel/23.1.0\n'
    txt = txt + 'source /scratch1/09778/jhodges/fds/FDS6VARS.sh\n'
    txt = txt + '\n'
    txt = txt + '# Launch\n'
    runningCount = 0
    for n, c in zip(numberOfProcesses, chids):
        txt = txt + 'ibrun -n %d -o %d task_affinity fds_impi_intel_linux %s.fds &\n'%(n, runningCount, c)
        runningCount = runningCount + n
    txt = txt + 'wait\n'
    return txt

if __name__ == "__main__":
    
    cases = [[1, 150, 1.1],
             [2, 200, 1.1],
             [3, 300, 1.1],
             [4, 400, 1.1],
             [5, 100, 0.6],
             [6, 200, 0.6],
             [7, 300, 0.6],
             [8, 400, 0.6],
             ]
    
    setnum = 0
    namespace = 'lat'
    (queue, runTime) = ('small', '24:00:00')
    for dx in [0.050, 0.025]:
        numberOfProcesses = []
        chids = []
        for case in cases:
            for convmodel in ['jet']:#['jet', 'default']:
                testnum = case[0]
                HRR = case[1]
                H = case[2]
                
                txt, numProcesses = buildInputFile(testnum, HRR, H, convmodel, dx)
                chid = getChid(testnum, convmodel, dx)
                outdir = "../"
                try:
                    os.mkdir(outdir)
                except:
                    pass
                with open(outdir + os.sep + chid+'.fds', 'w') as f:
                    f.write(txt)
                numberOfProcesses.append(numProcesses)
                chids.append(chid)
