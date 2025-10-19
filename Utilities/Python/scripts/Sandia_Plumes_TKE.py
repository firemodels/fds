
# Computes Turbulent Kinetic Energy (TKE) for Sandia Plumes Tests 24 and 35 and writes results to a csv file

import pandas as pd
import numpy as np
import os
import fdsplotlib

outdir = '../../../out/Sandia_Plumes/'

filename_out = [
    ['Sandia_CH4_1m_Test24_dx6cm_TKE.csv'  ,'Sandia_H2_1m_Test35_dx6cm_TKE.csv'],
    ['Sandia_CH4_1m_Test24_dx3cm_TKE.csv'  ,'Sandia_H2_1m_Test35_dx3cm_TKE.csv'],
    ['Sandia_CH4_1m_Test24_dx1p5cm_TKE.csv','Sandia_H2_1m_Test35_dx1p5cm_TKE.csv']
]

fds_line_file = [
    ['Sandia_CH4_1m_Test24_dx6cm_line.csv'  ,'Sandia_H2_1m_Test35_dx6cm_line.csv'],
    ['Sandia_CH4_1m_Test24_dx3cm_line.csv'  ,'Sandia_H2_1m_Test35_dx3cm_line.csv'],
    ['Sandia_CH4_1m_Test24_dx1p5cm_line.csv','Sandia_H2_1m_Test35_dx1p5cm_line.csv']
]


for j in range(3):  # resolution loop
    for k in range(2):  # test loop
        df = pd.read_csv(os.path.join(outdir,fds_line_file[j][k]),skiprows=1)
        x = df.iloc[:,0].values
        col = df.iloc[:,1:4].values
        for i in range(3):  # height loop
            w_rms = df.iloc[:,7+i].values
            u_rms = df.iloc[:,10+i].values
            col[:,i]= np.array(0.5*(w_rms**2+u_rms**2))
        A = pd.DataFrame({'x': x,
                          'TKE_p3': col[:,0],
                          'TKE_p5': col[:,1],
                          'TKE_p9': col[:,2]})
        A.to_csv(outdir + filename_out[j][k], index=False)


