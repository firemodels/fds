#6/2-21/2022 Noelle Crump, to convert data from Deep Fuel Bed Exp csv file
#into a paramfile to run through swaps.py to create input files for DFB
#code based on The MATLAB script of the same purpose by Ruddy Mell

import math
import pandas as pd
import os

    # read input data to dataframe
df = pd.read_csv('exp_params.csv', header=0)

    # initalize matrix to write final dataframe
FINAL = []

# Read data from the exp_params.csv file
for irow in df.index:
#    file_no = str(df.loc[irow,'RUN_NO'])
    burn_no = str(int(df.loc[irow,'BURN_NO']))    # burn number
    spacing = df.loc[irow,'SPACING']/100          # spacing, cm -> m
    spacing_s = str(int(df.loc[irow,'SPACING']))  # spacing as a string (without decimal)
    angle_d =  int(df.loc[irow,'SLOPE'])          # angle in degrees
    angle_s = str(angle_d)                        # angle as a string   (without decimal)
    angle_r = math.radians(df.loc[irow,'SLOPE'])    # angle in radians
    depth = df.loc[irow,'DEPTH']*0.0250           # in -> m !!!NOTE, in the exp the heights were conducted in cm,
                                                  # but are recorded in the .csv as in, ergo: *.0250 not *.0254
    depth_s = str(int(depth*100))                 # depth as a string in cm  (without decimal)
    m_ct = df.loc[irow,'MOISTURE_CONTENT']/100    # moisture, whole number percent to decimal percent
    res1 = df.loc[irow,'DX1']                     # res 1 to be used for platform spacing
    res2 = df.loc[irow,'DX2']                     # res 2 to be used for IJK
    s_res1 = 'p'+str(res1)[2:len(str(res1))]
    s_res2 = str(int(res2))
    temp = (df.loc[irow,'TMPA']-32.0)/1.8         # ambient temprature F -> C
    hmdy_rh = df.loc[irow,'HUMIDITY']             # ambient humidity
    burned_yn = int(df.loc[irow,'SPREAD'])        # did the fire spread 1=y 0=n.
    pct_burn = str(df.loc[irow,'BURN_PERCENT'])   # spread percentage (as whole %)
    param7IJK = res2                              # mesh IJK. later they should all be connected

# Record the baisic parameters:
    # begin the row for that file, with the filename:
    param1CHID = "burn"+burn_no+"_"+depth_s+"D_"+angle_s+"S_"+spacing_s+"L"
    param_line =[param1CHID+".fds"]
    # CHID and coments about the exp outcome
    param2TITLE = "USFS Deep Fuel Beds burn"+burn_no+" - "+depth_s+" cm deep - "+angle_s+" degree slope - "+spacing_s+" cm spacing"
    burnstring = ['did not burn','burned sucsessfully']
    burnpctstr = ['',' with a burn percentage of ' + pct_burn]
    param3BURN = "During the trial this run " + burnstring[burned_yn] + burnpctstr[burned_yn]
    param_line = param_line + [param1CHID] + [param2TITLE] + [param3BURN]
    # FDS params tempature, humidity, and moisture content
    param4TEMP = str(round(temp,2))
    param5HMDY = str(hmdy_rh)
    param6MRCT = str(round(m_ct,4))
    param_line = param_line + [param4TEMP] + [param5HMDY] + [param6MRCT]

    param_line = param_line + [str(param7IJK)] + [str(angle_d)]
    n_baisicparams = len(param_line)

# Record the fuel array parameters
    # calculate effective geometry descriptors
    startx = 0.6; # leftmost edge of the fuel rod arrays (and the platform they are on).
    rg_gap = spacing * math.cos(angle_r); # the slope distance, the effective horizontal distance between fuel rods.
    dz_ht = rg_gap * math.tan(angle_r)
#    eff_ht = rg_gap * math.tan(angle_r) + depth ; # effective height of bed
    rod_width = 0.15; # width of a fuel rod in m. (actually, should be .1524 but good enough for the cell size)
#    firegrid_top = max(round((max_ht+0.5),0)+2.0,2.0); # don't go any lower than 2.0 m.

    # calculate fuel array parameters
    xpos1 = startx
    xpos2 = xpos1 + rod_width
    if (angle_r == 0.0):
        zpos1 = 0
    else:
        zpos1 = (0.5*(xpos1+xpos2)-startx)*math.tan(angle_r)
    zpos2 = zpos1 + depth
    param_group = [rg_gap] + [(math.floor(4/spacing)-1)] + [zpos1] + [zpos2] + [dz_ht]
    param_line = param_line + param_group

    #round off large trailing decimals
    for n in range(n_baisicparams, len(param_line)):
        param_line[n] = round(param_line[n],4)

    #add paramline to FINAL matrix for each irow
    FINAL = FINAL + [param_line]

# Write paramfile.csv
    #make the header list
topline = ['Template.fds']
    #add basic params
for i in range(n_baisicparams-1):
    topline = topline + ["param"+str(i+1)]
    #add fuel params
for n in range(len(param_group)):
    topline = topline + ["paramf"+str(n+1)]
    #construct output dataframe
dfout = pd.DataFrame(FINAL, columns=topline)
    #write paramfile
dfout.to_csv('paramfile.csv', index=False)

# Build Input files, run swaps.py
os.system('python ../../../Utilities/Input_File_Tools/swaps.py')