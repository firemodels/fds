#!/usr/bin/env python

#generates 100 mesh inputs for Extinction.fds verification case
# K. Overholt, C. Weinschenk
# 9/2012


from __future__ import division

import numpy as np

fuel_seed= np.array([0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05])
temp_seed = np.array([300,475,650,825,1000,1175,1350,1525,1700,1875,300,475,650,825,1000,1175,1350,1525,1700,1875,300,475,650,825,1000,1175,1350,1525,1700,1875,300,475,650,825,1000,1175,1350,1525,1700,1875,300,475,650,825,1000,1175,1350,1525,1700,1875,300,475,650,825,1000,1175,1350,1525,1700,1875,300,475,650,825,1000,1175,1350,1525,1700,1875,300,475,650,825,1000,1175,1350,1525,1700,1875,300,475,650,825,1000,1175,1350,1525,1700,1875,300,475,650,825,1000,1175,1350,1525,1700,1875])

first_param = fuel_seed
second_param = 0.23 * (1 - fuel_seed)
temp = temp_seed-273.15

xb = []
for i in range(100):
    for j in range(1,21,2):
        z_dim = 0.1 * j
        for k in range(1,21,2):
            x_dim = 0.1 * k
            xb = np.append(xb,"%0.1f,%0.1f,0.0,0.1,%0.1f,%0.1f" % (x_dim-0.1, x_dim, z_dim-0.1, z_dim))

devc_numbers = range(1,101)

STRING="""
&MESH IJK=4,4,4, XB=%(xb)s /
&INIT XB=%(xb)s, SPEC_ID(1)='METHANE', MASS_FRACTION(1) = %(first_param)s, SPEC_ID(2)='OXYGEN', MASS_FRACTION(2)= %(second_param)s, TEMPERATURE = %(temp)s /
&ZONE XB=%(xb)s, /
&DEVC XB=%(xb)s, QUANTITY='HRRPUV', STATISTICS='VOLUME INTEGRAL', ID='HRR%(devc)s' /
&DEVC XB=%(xb)s, QUANTITY='TEMPERATURE', STATISTICS='MEAN', ID='TEMP%(devc)s' /
&DEVC XB=%(xb)s, QUANTITY='MASS FRACTION', STATISTICS='MEAN', SPEC_ID='OXYGEN', ID='O2_%(devc)s' /
&DEVC XB=%(xb)s, QUANTITY='MASS FRACTION', STATISTICS='MEAN', SPEC_ID='METHANE', ID='METHANE_%(devc)s' /"""

for i in range(100):
    print STRING % {'xb':xb[i], 'first_param':first_param[i], 'second_param':second_param[i], 'temp':temp[i], 'devc':devc_numbers[i]}

