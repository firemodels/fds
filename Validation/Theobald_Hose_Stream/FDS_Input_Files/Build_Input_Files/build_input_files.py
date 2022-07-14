# 7/13/2022 Noelle Crump. Converts data from theobald_effect_1981_fds.csv
# into a paramfile to run through swaps.py to build input files.
# Based on the build_input_files script of the same purpose in Deep_Fuel_Beds.

import math
import pandas as pd
import os

# read input data to dataframe
df = pd.read_csv('theobald_effect_1981_fds.csv', header=0)

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
    CHID = "Theobald_Test_"+file_num #+"_nozzle"+noz_name+"_"+rep_0(noz_D)+"D_"+rep_0(noz_P)+"P_"+noz_A+"A"
    filename = CHID + ".fds"
    T_str = "Theobald 1981 hose stream test series. nozzle"
    TITLE = T_str+noz_name+"- "+ noz_D+"mm diameter- "+noz_P+"bar pressure- "+noz_A+" degree angle"

    ORIENT3 = round(math.tan(math.radians(noz_angle)),2)           # z of orientation vector (x is 1)
    partvelo = 3.71*math.sqrt(noz_pres*14.504)              # velocity in m/s, (bar->psi before calc)
    flowrate = round((partvelo*math.pi*(noz_dmtr/2.0)**2)*60000.0,2) # flow rate in L/min (from m3/s)
    partvelo = round(partvelo,3)
    paramline = [filename]+[CHID]+[TITLE]+[partvelo]+[flowrate]+[ORIENT3]+[max_range]+[max_range+.1]

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

# move input files up one level
os.system('mv Theobald_Test_*.fds ../.')
