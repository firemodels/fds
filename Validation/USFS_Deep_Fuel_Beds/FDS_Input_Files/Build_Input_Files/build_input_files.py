#6/2-21/2022 Noelle Crump, to convert data from Deep Fuel Bed Exp csv file
#into a paramfile to run through swaps.py to create input files for DFB
#code based on The MATLAB script of the same purpose by Ruddy Mell

import math
import pandas as pd
import os

    # read input data to dataframe
df = pd.read_csv('../../../../../exp/USFS_Deep_Fuel_Beds/exp_params.csv', header=0)

    # initalize matrix to write final dataframe
FINAL = []

# Read data from the exp_params.csv file
for irow in df.index:
#    file_no = str(df.loc[irow,'RUN_NO'])
    burn_no = str(int(df.loc[irow,'BURN_NO']))    # burn number
    spacing = df.loc[irow,'SPACING']*0.0254       # spacing, in -> m
    spacing_s = str(int(spacing*98.4252))         # spacing as a string in cm (intervals 5cm)
    angle_d =  int(df.loc[irow,'SLOPE'])          # angle in degrees
    angle_s = str(angle_d)                        # angle as a string
    angle_r = math.radians(df.loc[irow,'SLOPE'])  # angle in radians
    depth = df.loc[irow,'DEPTH']*0.0254           # in -> m 
    depth_s = str(int(depth*98.4252))             # depth as a string in cm (121 -> 120)
    m_ct = df.loc[irow,'MOISTURE_CONTENT']/100    # moisture, whole number % to decimal %
    res1 = df.loc[irow,'DX1']                     # res 1 to be used for platform spacing
    res2 = df.loc[irow,'DX2']                     # res 2 to be used for IJK
    s_res1 = 'p'+str(res1)[2:len(str(res1))]      # res 1 as a string (0.25-> p25)
    s_res2 = str(int(res2))                       # res 2 as a string
    temp = (df.loc[irow,'TMPA']-32.0)/1.8         # ambient temprature F -> C
    hmdy_rh = df.loc[irow,'HUMIDITY']             # ambient humidity
    burned_yn = int(df.loc[irow,'SPREAD'])        # did the fire spread 1=y 0=n.
    pct_burn = int(df.loc[irow,'BURN_PERCENT'])   # spread percentage (as whole %)
    type_burn = int(round((pct_burn+49)/100,0)) + burned_yn # bad=0 marginal=1 good=2

# Record the baisic parameters:
    # begin the row for that file with the filename:
    param1CHID = "burn"+burn_no+"_"+depth_s+"D_"+angle_s+"S_"+spacing_s+"L"
    param_line =[param1CHID+".fds"]
    # CHID and coments about the exp outcome
    param2TITLE = "USFS Deep Fuel Beds burn"+burn_no+" - "+depth_s+" cm deep - "+angle_s+" degree slope - "+spacing_s+" cm spacing"
    burnstring = ['an unsuccessful burn','a marginal burn','a successful burn']
    burnpctstr = ['.',' with a burn percentage of ' + str(pct_burn)+'.','.']
    param3BURN = 'In EXP this run is clasified as ' + burnstring[type_burn] + burnpctstr[type_burn]
    param_line = param_line + [param1CHID] + [param2TITLE] + [param3BURN]
    # FDS params tempature, humidity, and moisture content
    param4TEMP = str(round(temp,2))
    param5HMDY = str(hmdy_rh)
    param6MRCT = str(round(m_ct,4))
    param_line = param_line + [param4TEMP] + [param5HMDY] + [param6MRCT]

    param_line = param_line + [str(res2)] + [str(angle_d)]
    n_basicparams = len(param_line)

# Record the fuel array parameters
    # calculate effective geometry descriptors
    startx = 0.6                            # leftmost edge of the fuel rod platform.
    rg_gap = spacing * math.cos(angle_r)    # horizontal distance between fuel rods.
    dz_ht = rg_gap * math.tan(angle_r)      # vertical distance between fueul rods
    rod_width = 0.1524                      # width of a fuel rod in m. (6 in)
    raise_part = 0.0762                     # bottom of rods about 3 in above fuel bed

    # calculate fuel array parameters
    xpos1 = startx
    xpos2 = xpos1 + rod_width
    if (angle_r == 0.0):
        zpos1 = 0 + raise_part
    else:
        zpos1 = (0.5*(xpos1+xpos2)-startx)*math.tan(angle_r) + raise_part
    #zpos2 = zpos1 + depth
    param_group = [rg_gap] + [dz_ht] + [(math.floor(4.8/spacing)-1)] + [zpos1] + [depth]

    # calculate obst array parameters
    block_base = 0.2                         # horizontal distance across each obst block.
    bed_base = 4.8 * math.cos(angle_r)       # horizontal distance underneath the tilted fuel bed
    dz_obst = block_base * math.tan(angle_r) # vertical distance between obst block bases
    param_group = param_group + [dz_obst] + [math.floor(bed_base/block_base) -1]

    param_line = param_line + param_group

    #round off large trailing decimals
    for n in range(n_basicparams, len(param_line)):
        param_line[n] = round(param_line[n],4)

    #add paramline to FINAL matrix for each irow
    FINAL = FINAL + [param_line]

# Write paramfile.csv
    #make the header list
topline = ['Template.fds']
    #add basic params
for i in range(n_basicparams-1):
    topline = topline + ["param"+str(i+1)]
    #add fuel params
for n in range(len(param_group)):
    topline = topline + ["paramf"+str(n+1)]
    #construct output dataframe
dfout = pd.DataFrame(FINAL, columns=topline)
    #write paramfile
dfout.to_csv('paramfile.csv', index=False)

# Build input files: run swaps.py
os.system('python ../../../../Utilities/Input_File_Tools/swaps.py')

# Move inpupt files up one level to FDS_Input_Files
os.system('mv ./burn*.fds ../.')