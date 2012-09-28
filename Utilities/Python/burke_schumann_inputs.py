#!/usr/bin/env python

#generates input file text for 15 mixture fraction reactors for Species/burke_schumann.fds verification case
# C. Weinschenk
# 9/2012

from __future__ import division

import numpy as np

fuel_seed= np.array([0,0.01438625,0.02836445,0.04195171,0.05518546,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])

first_param = fuel_seed

xyz = []
xb = []
lb = []
rb = []
fr = []
bk = []
tp = []

for j in range(1,31,2):
    x_dim = 0.1 * j
    xb = np.append(xb,"%0.1f,%0.1f,0.0,0.1,0.0,0.1" % (x_dim-0.1, x_dim))
devc_numbers = range(1,16)

STRING="""
&MESH IJK=4,4,4, XB=%(xb)s /
&ZONE XB=%(xb)s, /
&INIT XB=%(xb)s, SPEC_ID(1)='METHANE', MASS_FRACTION(1) = %(first_param)s /
&DEVC XB=%(xb)s, QUANTITY='TEMPERATURE', STATISTICS='MEAN', ID='TEMP_%(devc)s' /
&DEVC XB=%(xb)s, QUANTITY='DENSITY', STATISTICS='MEAN', ID='RHO_%(devc)s' /
&DEVC XB=%(xb)s, QUANTITY='SPECIFIC ENTHALPY', STATISTICS='MEAN', ID='H_%(devc)s' /
&DEVC XB=%(xb)s, QUANTITY='HRRPUV', STATISTICS='MEAN', ID='HRRPUV_%(devc)s' /
&DEVC XB=%(xb)s, QUANTITY='PRESSURE', STATISTICS='MEAN', ID='PRES_%(devc)s' /
&DEVC XB=%(xb)s, QUANTITY='MIXTURE FRACTION', STATISTICS='MEAN', ID='Z_%(devc)s' /
&DEVC XB=%(xb)s, QUANTITY='MASS FRACTION', STATISTICS='MEAN', SPEC_ID='OXYGEN', ID='O2_%(devc)s' /
&DEVC XB=%(xb)s, QUANTITY='MASS FRACTION', STATISTICS='MEAN', SPEC_ID='METHANE', ID='CH4_%(devc)s' /
&DEVC XB=%(xb)s, QUANTITY='MASS FRACTION', STATISTICS='MEAN', SPEC_ID='WATER VAPOR', ID='H2O_%(devc)s' /
&DEVC XB=%(xb)s, QUANTITY='MASS FRACTION', STATISTICS='MEAN', SPEC_ID='CARBON DIOXIDE', ID='CO2_%(devc)s' /
&DEVC XB=%(xb)s, QUANTITY='MASS FRACTION', STATISTICS='MEAN', SPEC_ID='NITROGEN', ID='N2_%(devc)s'  /"""

for i in range(15):
    print STRING % {'xb':xb[i], 'first_param':first_param[i],'devc':devc_numbers[i]}
