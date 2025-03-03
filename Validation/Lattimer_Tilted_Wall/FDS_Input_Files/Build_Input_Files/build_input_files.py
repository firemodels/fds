#6/1/2023 Jonathan Hodges, to generate titled wall input files based on the
#paramfile

import numpy as np

def generate_device_lines(theta, hght_list):
    lines = []     # generate devices corresponding to angle
    count = 0
    xs = []
    zs = []
    for height in hght_list:
        count += 1
        act_value = f"TC{count:02d}"
        hf_value = f"GHF{count:02d}"
        x_value = -height * np.sin(np.deg2rad(theta))
        z_value = height - (height * (1 - np.cos(np.deg2rad(theta))))
        line = f"&DEVC ID='{act_value}', QUANTITY='GAS TEMPERATURE', XYZ={x_value},0,{z_value}, IOR=1/\n&DEVC ID='{hf_value}', QUANTITY='GAUGE HEAT FLUX', XYZ={x_value},0,{z_value}, IOR=1/\n"
        lines.append(line)
        xs.append(x_value)
        zs.append(z_value)
    line="&DEVC XBP=0.0,%0.4f,0.0,0.0,0.0,%0.4f, QUANTITY='GAUGE HEAT FLUX', ID='GHF', POINTS=61, IOR=1, STATISTICS_START=50. /"%(np.nanmin(xs),np.nanmax(zs))
    lines.append(line)
    line="&DEVC XBP=0.0,%0.4f,0.0,0.0,0.0,%0.4f, QUANTITY='GAS TEMPERATURE', ID='TGAS', POINTS=61, IOR=1, STATISTICS_START=50., HIDE_COORDINATES=T /\n"%(np.nanmin(xs),np.nanmax(zs))
    lines.append(line)
    return lines

with open('paramfile.csv', 'r') as f:
    lines = f.readlines()

names  = [x.split(',')[0] for x in lines[1:]]
hrrs   = [float(x.split(',')[1]) for x in lines[1:]]
angles = [float(x.split(',')[2].replace('\n','')) for x in lines[1:]]

hght_list = [x / 100 for x in range(15, 140, 15)]  #board height

xs = [-x * np.sin(np.deg2rad(11)) for x in hght_list]
zs = [x - (x * (1 - np.cos(np.deg2rad(11)))) for x in hght_list]
hrrpuas = [x/(0.3*0.3) for x in hrrs]

templateFile = lines[0].split(',')[0]
with open(templateFile, 'r') as f:
    baseText = f.read()

for name, hrrpua, angle in zip(names, hrrpuas, angles):
    fileText = str(baseText)
    if '_obst' in name:
        walltxt = "&OBST ID='WALL', SURF_ID='ceramic board', XB=-0.0254,0.0,-0.30,0.30,0.0,1.5, /\n"
    else:
        walltxt = "&MOVE ID='SAMPLE-ROTATE', AXIS=0,-1,0, ROTATION_ANGLE, X0=0.0, Y0=0.0, Z0=0.0, DZ=0., DX=0., DY=0./\n"
        walltxt = walltxt + "&GEOM ID='WALL', SURF_ID='ceramic board', XB=-0.0254,0.0,-0.30,0.30,0.0,1.5, MOVE_ID='SAMPLE-ROTATE', CELL_BLOCK_IOR=-3,  /"
    fileText = fileText.replace('WALLSTRINGS',walltxt)
    
    if '_fine' in name:
        meshtxt = "&MESH ID='Mesh01', IJK=22,28,20, XB=-1.6,-1.325,-0.35,0.0,-0.2,0.05, MULT_ID='MESH'/\n"
        meshtxt = meshtxt + "&MULT ID='MESH', I_UPPER=7, J_UPPER=1, K_UPPER=7, DX=0.275, DY=0.35, DZ=0.25 /\n"
    elif '_coarse' in name:
        meshtxt = "&MESH ID='Mesh01', IJK=22,14,20, XB=-1.6,-0.5,-0.35,0.35,-0.2,0.8, MULT_ID='MESH'/\n"
        meshtxt = meshtxt + "&MULT ID='MESH', I_UPPER=1, K_UPPER=1, DX=1.1, DZ=1.0 /\n"
    else:
        meshtxt = "&MESH ID='Mesh01', IJK=44,28,40, XB=-1.6,-0.5,-0.35,0.35,-0.2,0.8, MULT_ID='MESH'/\n"
        meshtxt = meshtxt + "&MULT ID='MESH', I_UPPER=1, K_UPPER=1, DX=1.1, DZ=1.0 /\n"
    fileText = fileText.replace('MESHSTRINGS',meshtxt)
    
    fileText = fileText.replace('ROTATION_ANGLE','ROTATION_ANGLE=%0.1f'%(angle))
    fileText = fileText.replace('HRRPUA','HRRPUA=%0.1f'%(hrrpua))
    fileText = fileText.replace('CHID',"CHID='%s'"%(name.replace('.fds','')))
    device_lines = generate_device_lines(angle, hght_list)
    fileText = fileText + '\n'.join(device_lines)
    fileText + '\n&TAIL /\n'
    
    
        
    with open("..//%s"%(name), 'w') as f:
        f.write(fileText)
