#!/usr/bin/python
#McDermott
#08 Sep 2017
#
# This script is handy for figuring out the base mesh in a MULT arrangement.
#
# Inputs:
#   IJK   = single mesh equivalent IJK
#   XB    = single mesh equivalent XB
#   MBLKS = mesh block arrangment
#   M0    = indices of the "base" block for MULT

import math

# input FDS mesh line information

IJK=221,221,241
XB=-1.5,1.5,-1.5,1.5,-.3,3
MBLKS=5,5,5
M0=3,3,1

DX=(XB[1]-XB[0])/IJK[0]
DY=(XB[3]-XB[2])/IJK[1]
DZ=(XB[5]-XB[4])/IJK[2]

# write out mesh block M0

NX = math.ceil(IJK[0]/MBLKS[0])
NY = math.ceil(IJK[1]/MBLKS[1])
NZ = math.ceil(IJK[2]/MBLKS[2])
LX = NX*DX
LY = NY*DY
LZ = NZ*DZ

print(DX,DY,DZ)
print(LX,LY,LZ)

XB0 = (-0.5*LX, 0.5*LX, -0.5*LY, 0.5*LY, XB[4], XB[4]+LZ)

IJK_LOC=str((NX,NY,NZ))
XB_LOC=str(XB0)

print('&MESH IJK='+IJK_LOC[1:-1]+', XB='+XB_LOC[1:-1]+'/')