#6/2/2022 Noelle Crump, to convert data from Deep Fuel Bed Exp csv file
#into a paramfile to run through swaps.py to create input files for DFB
#code based on The MATLAB script of the same purpose by Rudy

import math
import csv
import os
import numpy

file_n = 'experiments.csv'
paramfile = 'paramfile.csv'

n_data = numpy.loadtxt(fname= file_n, delimiter=',')
data = n_data
#initailize final matrix to write excluding the top line
FINAL = []

#Initalize data to write param file:    
        #read the data in the file:
for irow in range(len(data)):    
    file_no = str(data[irow,0])             #
    burn_no = str(int(data[irow,1]))        # burn number
    spacing = data[irow,2]/100              # spacing, cm -> m
    spacing_s = str(int(data[irow,2]))      # spacing as a string (without decimal)
    angle = data[irow,3]
    angle_s = str(int(data[irow,3]))       # angle as a string   (without decimal)
    angle_r = numpy.radians(data[irow,3])  # degrees-> radians
    depth = data[irow,4]*0.0254            # in -> meters
    depth_s = str(int(data[irow,4]))       # depth as a string   (without decimal)
    m_ct = data[irow,5]/100                # percent to decimal percent                                     
    res1 = data[irow,6]                    # res 1 to be used for platform spacing
    res2 = data[irow,7]                    # res 1 to be used for IJK
    temp = (data[irow,8]-32)/1.8           # ambient temprature F -> C
    hmdy_rh = data[irow,9]                 # ambient humidity
    burned_yn =int(data[irow,10])          # did the fire spread 1=y 0=n.
    pct_burn =str(data[irow,11])+'%.'      # spread percentage (as whole %)
    #roundint = int()                      #for readability all data is rounded before it prints to the param file
                                           #the minimum round number is 4, but at higer resolutions, it should increase
        
        #Hardwires:
    res1= 0.25          # temporary hardwire for low resoluiton runs
    n_ptcls = 1000      # for now, hardwire Number of particles and
    param7IJK = res2      # mesh IJK. later they should all be connected


#Start Writing the file: 
        #Begin the row for that file, with the Filename:
    param_line =[burn_no+"burnno_"+depth_s+"D_"+angle_s+"S_"+spacing_s+"L_"+"DFB.fds"] 
        # HEAD, CHID, and coments about the exp outcome
    param1HEAD =depth_s+"D_"+angle_s+"S_"+spacing_s+"L_"+"res-"+s_res
    param2TITLE ="'Deep fuel bed-"+depth_s+" in deep- "+angle_s+" degrees- "+spacing_s+" cm spacing. expburn#"+burn_no+". Matlab RUN"+file_no+"'"
    burnstring=['did not burn','burned sucsessfully']
    burnpctstr =['.',' with a burn percentage of '+pct_burn ]
    param3BURN = "During the trial this run "+ burnstring[burned_yn]+burnpctstr[burned_yn]
    param_line = param_line + [param1HEAD]+[param2TITLE]+[param3BURN]
        # FDS Params tempature, humidity, and moisture content
    param4TEMP = str(round(temp,2))
    param5HMDY = str(hmdy_rh)
    param6MRCT = str(round(m_ct,4))
    param_line = param_line +[param4TEMP]+[param5HMDY]+[param6MRCT]
    
    #~-Fuel particle density =  514 kg/m^3 from the notes(Russ/Mark)
   # n_ptcls = str(1000)# for now, hardwire Number of particles
    #param7IJK = 30 #and mesh resolution. later they should be connected
    param_line = param_line+[str(param7IJK)]
    n_baisicparams= len(param_line)
    
#Writing the arrays:
        #The fuel arrays:
    # calculate effective geometry descriptors
    startx = 0.6; # leftmost edge of the fuel rod arrays (and the platform they are on).
    rg_gap = spacing * math.cos(angle_r); # the slope distance, the effective horizontal distance between fuel rods.
    eff_ht = rg_gap * math.tan(angle_r) + depth ; # effective height of bed
    rod_width = 0.15; # width of a fuel rod in m. (actually, should be .1524 but good enough for the cell size)
    max_ht = 4.8*math.tan(angle_r)+depth ; # this is the highest point that fuel goes to.
    firegrid_top = max(round((max_ht+0.5),0)+2.0,2.0); # don't go any lower than 2.0 m.
    
    #make the x vect
    endx = startx + 4.8 - 0.15
    s_xvect = math.ceil((endx-startx)/spacing)
    xvect=s_xvect*[startx]
    for i in range(s_xvect):
        xvect[i]=round(xvect[i]+(i*spacing),4)
        
    # make the fuel bed coordinates
    param_group = []
    for ind in range(2):
        xpos1 = xvect[(ind)]
        xpos2 =xpos1 +rod_width
        zpos1 = (0.5*(xpos1+xpos2)-xvect[0])*math.tan(angle_r) -0.3

        if (angle_r == 0.0):
            zpos2 = depth
            zpos1 = 0.0
        else:
            zpos1 = zpos1 + 0.3
            zpos2 = zpos1 + depth
 
        paramf1 =xpos1
        paramf2 =xpos2
        paramf3 =zpos1
        paramf4 =zpos2
        paramf5 = n_ptcls
        if(ind==0):
            param_group = [paramf1]+[paramf2]+[paramf3]+[paramf4]+[paramf5]
            vent_end = xpos2 #save for vent for later
        elif(ind==1):
            param_group = param_group + [paramf1 -param_group[0]] + [paramf3 -param_group[2]]+[s_xvect-1]
    param_line = param_line+param_group  

        #the OBST array:
    #make the o vect
    endx = startx + 4.8
    s_Oxvect = math.ceil((endx-startx)/res1)

    Oxvect=s_Oxvect*[startx]
    for i in range(s_Oxvect):
        Oxvect[i]=Oxvect[i]+(i*res1)
        
    #make the obst coordinates
    param_group = []
    for ind in range(2):
        Oxpos1 = Oxvect[ind]
        Oxpos2 = Oxpos1 +res1
        Ozpos1b = (0 - res1)
        Ozpos2 = ((Oxpos1-Oxvect[0])*math.tan(angle_r))
    
        if(ind==0):         
            paramo1 = Ozpos2 + .05*res1 #saving a vent height param
            param_group = [Oxpos1]+[Oxpos2]+[Ozpos1b]+[Ozpos2]
        elif(ind==1):
            param_group = param_group + [Oxpos1 -param_group[0]] + [Ozpos2 -param_group[3]] +[s_Oxvect-1]

    if(param_group[4]>.20):
        n_Vmult = 0
    else:
        n_Vmult =int(round((rod_width/param_group[4]),0)) -1
    param_group[0] = paramo1
    param_line = param_line+param_group +[n_Vmult]
    #replace paramo1 with n for vent
    
    
    for n in range(n_baisicparams, len(param_line)):
        param_line[n]=round(param_line[n],4)


        
    FINAL = FINAL + [param_line]
    
# Write the paramfile    
if os.path.isfile(paramfile):
    os.remove(paramfile)
print('writing to param file')

topline = ['Template.fds']

#add baisic params
for i in range((n_baisicparams-1)):
    topline = topline +["'param"+ str(i+1)+"'"]

#add tree params
param_header=[]
for n in range(8):
    topline = topline+ ["'paramf"+str(n+1)+"'"]

for m in range(8):
    topline = topline+ ["'paramo"+str(m+1)+"'"]

#write final param file
with open(paramfile, 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(topline)
    for row in range(len(data)):
        writer.writerow(FINAL[row])
