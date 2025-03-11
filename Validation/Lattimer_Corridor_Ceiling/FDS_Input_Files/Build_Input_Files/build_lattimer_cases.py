# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 18:21:38 2024

@author: jhodges
"""
import os

def getChid(testnum, dx):
    suffix1 = '_%03dmm'%(dx*1e3)
    chid = 'lattimer_test%0d%s'%(testnum, suffix1)
    return chid

def buildInputFile(testnum, HRR, H, dx):
    chid = getChid(testnum, dx)
    suffix2 = "HEAT_TRANSFER_MODEL='IMPINGING JET', CONVECTION_LENGTH_SCALE=0.77, "
        
    txt = "&HEAD CHID='%s', /\n\n"%(chid)
    txt = txt + "&TIME T_END=200., TIME_SHRINK_FACTOR=10 /\n"
    txt = txt + "&DUMP DT_HRR=10.,  /\n"
    
    if dx == 0.05: 
        numProcesses = 7
        txt = txt + "&MESH XB=-3.2,-2.7,-0.875,0.925,0.0,2.4, IJK=10,36,48, MULT_ID='MESH' /\n"
        txt = txt + "&MULT ID='MESH', DX=0.5, I_UPPER=6, /\n\n"
        npoints = 47
    if dx == 0.025:
        txt = txt + "&MESH XB=-3.2,-2.325,-0.8875,0.0125,0.0,1.2, IJK=35,36,48, MULT_ID='MESH' /\n"
        txt = txt + "&MULT ID='MESH', DX=0.875, DY=0.9, DZ=1.2, I_UPPER=3, J_UPPER=1, K_UPPER=1 /\n\n"
        numProcesses = 16
        npoints = 47
    if dx == 0.0125:
        txt = txt + "&MESH XB=-2.95,-2.2,-0.65,-0.325,0.0,0.5625, IJK=60,26,45, MULT_ID='MESH' /\n"
        txt = txt + "&MULT ID='MESH', DX=0.75, DY=0.325, DZ=0.5625, I_UPPER=3, J_UPPER=3, K_UPPER=3 /\n\n"
        numProcesses = 64
        npoints = 47
    dz = dx
    
    txt = txt + "&REAC FUEL='PROPANE', SOOT_YIELD=0.024, CO_YIELD=0.005 / \n\n"
    txt = txt + "&RADI C_MIN=1. /\n\n"
    
    txt = txt + "&MATL ID= 'CONCRETE', FYI = 'Drysdale, Intro to Fire Dynamics', SPECIFIC_HEAT = 0.90, DENSITY = 2000., CONDUCTIVITY  = 1.0, /\n"
    txt = txt + "&SURF ID = 'FLOOR', MATL_ID = 'CONCRETE', THICKNESS  = 0.1, COLOR = 'GRAY' /\n\n"
    
    txt = txt + "&MATL ID='CERAMIC-FIBER', SPECIFIC_HEAT=1.14, DENSITY=96, CONDUCTIVITY=0.14 /\n"
    txt = txt + "&SURF ID='CEILING', HEAT_TRANSFER_MODEL='IMPINGING JET', EMISSIVITY=0.8, MATL_ID='CERAMIC-FIBER', THICKNESS=0.0254, COLOR='GRAY 80' /\n"
    txt = txt + "&SURF ID='WALL', EMISSIVITY=0.8, MATL_ID='CERAMIC-FIBER', THICKNESS=0.0254, COLOR='GRAY 40' /\n\n"
    
    txt = txt + "&SURF ID='fire', HRRPUA=%0.4f, TMP_FRONT=400, COLOR='RED', /\n\n"%(HRR/(0.495))
    
    txt = txt + "&OBST ID='CEILING', XB=-2.4,%0.4f,%0.4f,%0.4f,2.1,%0.4f, SURF_IDS='WALL','WALL','CEILING', /\n"%(dx,-0.6-dx,0.6+dx,2.1+dz)
    txt = txt + "&OBST ID='BACK-WALL', XB=0.0,%0.4f,%0.4f,%0.4f,0.9,2.1, SURF_ID='WALL', /\n"%(0.0+dx, -0.6-dx, 0.6+dx)
    txt = txt + "&OBST ID='SIDEWALL-1', XB=-2.4,0.0,%0.4f,-0.6,1.5,2.1, SURF_ID='WALL', /\n"%(-0.6-dx)
    txt = txt + "&OBST ID='SIDEWALL-2', XB=-2.4,0.0,0.6,%0.4f,1.5,2.1, SURF_ID='WALL', /\n"%(0.6+dx)
    txt = txt + "&OBST ID='SIDEWALL-3', XB=-1.2,0.0,%0.4f,-0.6,0.9,1.5, SURF_ID='WALL', /\n"%(-0.6-dx)
    txt = txt + "&OBST ID='SIDEWALL-2', XB=-1.2,0.0,0.6,%0.4f,0.9,1.5, SURF_ID='WALL', /\n\n"%(0.6+dx)
    
    txt = txt + "&OBST ID='FIRE', XB=-0.45,0.0,-0.55,0.55,%0.4f,%0.4f, SURF_IDS='fire','WALL','WALL' /\n"%(2.1-H-0.15,2.1-H)
    txt = txt + "&OBST ID='BLOCK', XB=-0.45,0.0,-0.55,0.55,0.0,%0.4f, SURF_ID='FLOOR' /\n\n"%(2.1-H-0.15)
    
    txt = txt + "&VENT DB='XMIN', SURF_ID='OPEN' /\n"
    txt = txt + "&VENT DB='XMAX', SURF_ID='OPEN' /\n"
    txt = txt + "&VENT DB='YMIN', SURF_ID='OPEN' /\n"
    txt = txt + "&VENT DB='YMAX', SURF_ID='OPEN' /\n"
    txt = txt + "&VENT DB='ZMIN', SURF_ID='FLOOR', /\n"
    txt = txt + "&VENT DB='ZMAX', SURF_ID='OPEN' /\n\n"
    
    hide_coordinates = False
    suffix = ""
    for ID, QTY in zip(['GHF_profile','HTC_profile','CHF_profile'], ['GAUGE HEAT FLUX', 'HEAT TRANSFER COEFFICIENT', 'CONVECTIVE HEAT FLUX']):
        if hide_coordinates: suffix = ", HIDE_COORDINATES=T"
        txt = txt + "&DEVC ID='%s', XB=-2.35,-0.05,0.0,0.0,2.1,2.1, IOR=-3, POINTS=%d, QUANTITY='%s', STATISTICS_START=50.%s /\n"%(ID, npoints, QTY, suffix)
        hide_coordinates=True
    txt = txt + "\n"
    
    for qty in ['HEAT TRANSFER COEFFICIENT', 'GAUGE HEAT FLUX', 'CONVECTIVE HEAT FLUX', 'WALL TEMPERATURE']:
        txt = txt + "&BNDF QUANTITY='%s', "%(qty)
        txt = txt + "CELL_CENTERED=T /\n"
    
    for qty in ['TEMPERATURE']:
        txt = txt + "&SLCF PBY=0, QUANTITY='%s', "%(qty)
        txt = txt + "CELL_CENTERED=T /\n\n"
    txt = txt + '&TAIL /\n'
    
    return txt, numProcesses

if __name__ == "__main__":
    cases = [[1, 150, 1.1],
             [2, 200, 1.1],
             [3, 300, 1.1],
             [4, 400, 1.1],
             [5, 100, 0.6],
             [6, 200, 0.6],
             [7, 300, 0.6],
             [8, 400, 0.6]]
    outdir = ".." + os.sep
    try:
        os.mkdir(outdir)
    except:
        pass
    for dx in [0.050, 0.025, 0.0125]:
        for case in cases:
            testnum, HRR, H = case[0], case[1], case[2]
            txt, numProcesses = buildInputFile(testnum, HRR, H, dx)
            chid = getChid(testnum, dx)
            with open(outdir + os.sep + chid+'.fds', 'w') as f:
                f.write(txt)
