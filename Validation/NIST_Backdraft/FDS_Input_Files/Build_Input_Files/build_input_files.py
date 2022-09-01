# Based on the build_input_files script of the same purpose in Deep_Fuel_Beds.

import math
import pandas as pd
import os
import numpy as np

# read input data to dataframe
df = pd.read_csv('../../../../../exp/NIST_Backdraft/Backdraft_Two_Zone_init.csv', header=0)

# initalize matrix to write final dataframe
FINAL = []

# function gets rid of .0 s or replaces . with p
def rep_0(str_in):
    if str_in[2:]=='.0':
        str_out = str_in[:2]
    else:
        str_out = str_in.replace('.','p')
    return(str_out)

Grid_Size_Vect   = ['coarse','medium','fine']
Ignitor_Pos_Vect = ['low','mid']
Grid_Size        = "coarse"
Ignitor_Pos      = "low"
Grid_Size_File   = Grid_Size + "_mesh.fds"
Ignitor_Pos_File = "ignitor_" + Ignitor_Pos + ".fds"

# read data from the Backdraft_Two_Zone_init.csv file
for Grid_Size in Grid_Size_Vect:
    Grid_Size_File   = Grid_Size + "_mesh.fds"
    for Ignitor_Pos in Ignitor_Pos_Vect:
        Ignitor_Pos_File = "ignitor_" + Ignitor_Pos + ".fds"
        for irow in df.index:

            Fuel        = str(df.loc[irow,'Fuel'])
            Reac_File   = Fuel.lower() + "_reac.fds"
            Fire_Size   = str(df.loc[irow,'Fire_Size_kW'])
            Fire_Size   = rep_0(Fire_Size)
            Flow_Time   = str(df.loc[irow,'Fuel_Flow_Time_s'])
            Up_Temp     = str(df.loc[irow,'Upper_Temperature_C'])
            Up_X_Fuel   = str(df.loc[irow,'Upper_X_Fuel'])
            Up_X_O2     = str(df.loc[irow,'Upper_X_O2'])
            Up_X_CO2    = str(df.loc[irow,'Upper_X_CO2'])
            Up_X_CO     = str(df.loc[irow,'Upper_X_CO'])
            Up_X_H2O    = str(df.loc[irow,'Upper_X_H2O'])
            Up_X_N2     = str(df.loc[irow,'Upper_X_H2O'])
            Lo_Temp     = str(df.loc[irow,'Lower_Temperature_C'])
            Lo_X_Fuel   = str(df.loc[irow,'Lower_X_Fuel'])
            Lo_X_O2     = str(df.loc[irow,'Lower_X_O2'])
            Lo_X_CO2    = str(df.loc[irow,'Lower_X_CO2'])
            Lo_X_CO     = str(df.loc[irow,'Lower_X_CO'])
            Lo_X_H2O    = str(df.loc[irow,'Lower_X_H2O'])
            Lo_X_N2     = str(df.loc[irow,'Lower_X_H2O'])

            # make the row for that file
            CHID = "NIST_Backdraft_" + Fuel + "_" + Fire_Size + "kW_" + Flow_Time + "s_" + Ignitor_Pos + '_' + Grid_Size
            filename = CHID + ".fds"
            T_str = "NIST Backdraft experiment : "
            TITLE = "NIST Backdraft experiment : " + Fuel + "-" + Fire_Size + " kW-" + Flow_Time + "s-Ignitor position " + Ignitor_Pos + "."
            paramline = [filename] + [CHID] + [TITLE] + [Grid_Size_File] + [Ignitor_Pos_File] + [Reac_File] + [Fuel.upper()]
            paramline = paramline  + [Up_Temp] + [Up_X_Fuel] + [Up_X_O2] + [Up_X_CO2] + [Up_X_CO] + [Up_X_H2O] + [Up_X_N2]
            paramline = paramline  + [Lo_Temp] + [Lo_X_Fuel] + [Lo_X_O2] + [Lo_X_CO2] + [Lo_X_CO] + [Lo_X_H2O] + [Lo_X_N2]

            # add paramline to FINAL matrix for each irow
            FINAL = FINAL + [paramline]


# make the header list
topline = ["NIST_BackdraftTemplate.fds"] + ["paramCHID"] + ["paramTITLE"]
for i in range(3,21):
    topline = topline + ["param"+str(i)]

# make fdout and wirte paramfile
dfout = pd.DataFrame(FINAL, columns=topline)
dfout.to_csv('paramfile.csv', index=False)

# build input files: run swaps.py
os.system('python ../../../../Utilities/Input_File_Tools/swaps.py')

# Move inpupt files up one level to FDS_Input_Files
os.system('mv NIST_Backdraft_*.fds ../')
