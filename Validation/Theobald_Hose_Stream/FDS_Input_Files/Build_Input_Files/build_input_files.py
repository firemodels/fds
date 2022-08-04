# 7/13/2022 Noelle Crump. Converts data from theobald_effect_1981_fds.csv
# into a paramfile to run through swaps.py to build input files.
# Based on the build_input_files script of the same purpose in Deep_Fuel_Beds.

import math
import pandas as pd
import os

# read input data to dataframe
df = pd.read_csv('../../../../../exp/Theobald_Hose_Stream/theobald_effect_1981_fds.csv', header=0)

# initalize matrix to write final dataframe
FINAL = []

# function gets rid of .0 s or replaces . with p
def rep_0(str_in):
    if str_in[2:]=='.0':
        str_out = str_in[:2]
    else:
        str_out = str_in.replace('.','p')
    return(str_out)

# read data from the theobald_effect_1981_fds.csv file
for irow in df.index:
    noz_name = str(df.loc[irow,'nozzle']).replace(' ','')
    file_num = str(df.loc[irow,'Unnamed: 0'])
    noz_dmtr = (df.loc[irow,'nozzle diameter (mm)']/1000.0)   #(m)
    noz_pres = df.loc[irow,'dP (bar)']                      #(bar)
    noz_angle = df.loc[irow,'firing angle (deg.)']      #(degrees)
    max_range = df.loc[irow,'max range (m)']                  #(m)
    noz_D = str(df.loc[irow,'nozzle diameter (mm)'])
    noz_P = str(df.loc[irow,'dP (bar)'])
    noz_A = str(int(df.loc[irow,'firing angle (deg.)']))

# make the row for that file
    CHID = "Theobald_Test_" + file_num
    filename = CHID + ".fds"
    T_str = "Theobald 1981 hose stream test series. nozzle"
    TITLE = T_str+noz_name+"- "+noz_D+"mm diameter- "+noz_P+"bar pressure- "+noz_A+" degree angle"

    orient_z = round(math.tan(math.radians(noz_angle)),2)           # z of orientation vector (1,0,z)
    efx_velo = 3.71*math.sqrt(noz_pres*14.504)       # efflux velocity in m/s, (bar->psi before calc)
    flowrate = round((efx_velo*math.pi*(noz_dmtr/2.0)**2)*60000.0,2) # flow rate in L/min (from m3/s)
    We = 13864.0*noz_dmtr*efx_velo**2              # Weber Number (@ 20C, diameter(m), velocity(m/s))
    Re = 996605.43*noz_dmtr*efx_velo            # Reynolds Number (@ 20C, diameter(m), velocity(m/s))
    breaklength50 = 1.7*noz_dmtr*(We**0.5)*(Re*1e-4)**(-0.625) # length @ stream 50% discontinuous(m)

    fdsbreaklength = breaklength50*1.0        # FDS PRIMARY_BREAKUP_LENGTH (m)

    paramline = [filename] + [CHID] + [TITLE] + [round(efx_velo,3)] + [flowrate] + [orient_z]
    paramline = paramline  + [max_range] + [round((max_range+0.1),1)] + [round(fdsbreaklength,2)]

# add paramline to FINAL matrix for each irow
    FINAL = FINAL + [paramline]

# make the header list
topline = ['theobald_Template.fds']
for i in range(1, len(paramline)):
    topline = topline + ["param"+str(i)]

# make fdout and wirte paramfile
dfout = pd.DataFrame(FINAL, columns=topline)
dfout.to_csv('paramfile.csv', index=False)

# build input files: run swaps.py
os.system('python ../../../../Utilities/Input_File_Tools/swaps.py')

# Move inpupt files up one level to FDS_Input_Files
os.system('mv Theobald_Test_*.fds ../.')