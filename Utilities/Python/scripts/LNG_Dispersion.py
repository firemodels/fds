#!/usr/bin/env python
"""
McGrattan
12-6-2016
LNG_Disperion.m

Reads and transposes line data at various times in simulation
"""

import os
import numpy as np
import pandas as pd

outdir = '../../../out/LNG_Dispersion/'
expdir = '../../../exp/LNG_Dispersion/'

labels = ['Burro3','Burro7','Burro8','Burro9',
          'Coyote3','Coyote5','Coyote6',
          'Falcon1','Falcon3','Falcon4',
          'MaplinSands27','MaplinSands34','MaplinSands35']

c_start = [[2,20,32,41],[2,20,41,50],[2,20,41,62],[2,20,32],
           [2,29,44,50,65],[2,29,44,53,68],[2,29,44,53,68],
           [2,25,49],[2,21,44],[2,21,53],
           [2,5,8,12,13,19,22,25],[2,5],[2,5,8]]

c_end   = [[19,31,40,46],[19,40,49,52],[19,40,61,73],[19,31,40],
           [28,43,49,64,67],[28,43,52,67,73],[28,43,52,67,73],
           [24,48,65],[20,43,61],[20,52,67],
           [4,7,11,12,18,21,24,27],[4,7],[4,7,10]]

t_start = [[30,30,60,170],[40,40,80,140],[60,150,400,650],[40,80,140],
           [50,50,50,50,50],[40,40,40,40,40],[38,38,38,38,38],
           [100,125,150],[100,125,175],[100,150,150],
           [0,0,0,0,0,0,0,0],[0,0],[0,0,0]]

t_end   = [[130,130,160,270],[180,180,220,280],[140,230,480,730],[180,220,280],
           [100,100,100,100,100],[130,130,130,130,130],[108,108,108,108,108],
           [200,225,250],[250,275,325],[350,400,400],
           [160,160,160,160,160,160,160,160],[95,95],[135,135,135]]

H1 = ['m','%']
H2 = ['x','X_CH4']

for j in range(13):
    # Read devc and exp data (skip header lines to match MATLAB importdata)
    M = pd.read_csv(os.path.join(outdir, f"{labels[j]}_devc.csv"), skiprows=2)
    E = pd.read_csv(os.path.join(expdir, f"{labels[j]}_exp.csv"), skiprows=2)
    
    # Prepare output file
    with open(os.path.join(outdir, f"{labels[j]}.csv"), 'w') as fid:
        fid.write(f"{H1[0]}, {H1[1]}\n")
        fid.write(f"{H2[0]}, {H2[1]}\n")
        
        for i in range(len(E)):
            time = M.iloc[:,0].values
            idx = np.where((time > t_start[j][i]) & (time < t_end[j][i]))[0]
            if len(idx) == 0:
                val = 0.01
            else:
                data_segment = M.iloc[idx, c_start[j][i]-1:c_end[j][i]].values  # adjust index (MATLAB 1-based)
                val = max(0.01, 100. * np.max(data_segment))
            fid.write(f"{E.iloc[i,0]:5.1f}, {val:6.2f}\n")
