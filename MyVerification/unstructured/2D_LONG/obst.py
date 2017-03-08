#me -*- coding: utf-8 -*-

#####################################################################################
# Initial import
#####################################################################################
import sys
from pylab import *
from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib.font_manager as fnt



nmeshes = int(sys.argv[1])


xs = -0.8
xf =  0.8
ys = -0.01
yf =  0.01
zs = -0.8
zf =  0.8

for i in range(nmeshes):

   print "&OBST XB = %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f /" %(xs,xf,ys,yf,zs,zf)
   xs += 4.8
   xf += 4.8




