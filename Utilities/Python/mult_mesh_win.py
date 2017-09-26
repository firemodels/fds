#!/usr/bin/python

from __future__ import division
from decimal import *
from math import sqrt
import string

import math

def start() :
    print("\n------------------------------------------------------------------")
    print   "\t*** Welcome To Multy Mesh generator *** "
    print """\nAuthor: Salah Benkorichi\n26 Sep 2017\n\nThis is a handy script for figuring out the base mesh in a MULT arrangement.

    Inputs:
	  IJK   = single mesh equivalent IJK, exp: 221 221 241
	  XB    = single mesh equivalent XB,  exp: -1.5 1.5 -1.5 -3 3
	  MBLKS = mesh block arrangment,      exp: 5 5 5
	  M0    = indices of the "base" block for MULT, exp: 3 3 1
      \nDX,DY,DZ and mesh line will be output in mesh.txt file"""

    print("------------------------------------------------------------------")		

	#his lines ask for FDS mesh line information

    I,J,K =(raw_input("\nEnter IJK as I J K:  ").split())
    I=int(I)
    J=int(J)
    K=int(K)
    x1,x2,y1,y2,z1,z2 = (raw_input("\nEnter XB as x1 x2 y1 y2 z1 z2:  ").split())
    x1=float(x1)
    x2=float(x2)
    y1=float(y1)
    y2=float(y2)
    z1=float(z1)
    z2=float(z2)
    MBLKS_i,MBLKS_j,MBLKS_k = (raw_input("\nEnter MBLKS as MBLKS_i MBLKS_j MBLKS_k: ").split())
    MBLKS_i=int(MBLKS_i)
    MBLKS_j=int(MBLKS_j)
    MBLKS_k=int(MBLKS_k)
    M0_i,M0_j,M0_k = (raw_input("\nEnter M0 as M0_i M0_j M0_k: ").split())
    M0_i=int(M0_i)
    M0_j=int(M0_j)
    M0_k=int(M0_k)

    IJK=[I,J,K]
    XB=[x1,x2,y1,y2,z1,z2]
    MBLKS=[MBLKS_i,MBLKS_j,MBLKS_k]
    M0 = [M0_i,M0_j,M0_k]

    # write out mesh block M0

    DX=float((XB[1]-XB[0])/IJK[0])
    DY=float((XB[3]-XB[2])/IJK[1])
    DZ=float((XB[5]-XB[4])/IJK[2])

	# # write out mesh block M0

    NX = math.ceil(IJK[0]/MBLKS[0])
    NY = math.ceil(IJK[1]/MBLKS[1])
    NZ = math.ceil(IJK[2]/MBLKS[2])
    LX = NX*DX
    LY = NY*DY
    LZ = NZ*DZ

    print "DX= %s, DY= %s, DZ= %s" % (DX,DY,DZ)
    print "LX= %s, LY= %s, LZ= %s" % (LX,LY,LZ)

    XB0 = (-0.5*LX, 0.5*LX, -0.5*LY, 0.5*LY, XB[4], XB[4]+LZ)

    IJK_LOC=str((NX,NY,NZ))
    XB_LOC=str(XB0)

    mesh = '&MESH IJK='+IJK_LOC[1:-1]+', XB='+XB_LOC[1:-1]+'/'
    print mesh
    text_file = open("mesh.txt", "w")
    text_file.write("DX= %s, DY= %s, DZ= %s \nmesh: %s" % (DX,DY,DZ,mesh))
    text_file.close()
    return try_again()

    # Exit loop at the end of generation
def try_again():
    while True:
	
       exit_inp = raw_input("\nDo you want to try again? Y/N ")
       if exit_inp.lower() == ("n"):
          print "\nclosing!"
          import os
          quit(1)
       elif exit_inp.lower() == ("y"):
          return start()
       else:
          print "\nplease enter valid input"
                          
start()
