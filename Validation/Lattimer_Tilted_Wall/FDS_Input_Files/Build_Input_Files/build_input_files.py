#6/1/2023 Jonathan Hodges, to generate titled wall input files based on the
#paramfile

import numpy as np

def generate_device_lines(theta, hght_list):
    lines = []     # generate devices corresponding to angle
    count = 0
    for height in hght_list:
        count += 1
        act_value = f"TC{count:02d}"
        hf_value = f"GHF{count:02d}"
        x_value = -height * np.sin(np.deg2rad(theta))
        z_value = height - (height * (1 - np.cos(np.deg2rad(theta))))
        line = f"&DEVC ID='{act_value}', QUANTITY='GAS TEMPERATURE', XYZ={x_value},0,{z_value}, IOR=1/\n&DEVC ID='{hf_value}', QUANTITY='GAUGE HEAT FLUX', XYZ={x_value},0,{z_value}, IOR=1/\n"
        lines.append(line)
    return lines

with open('paramfile.csv', 'r') as f:
    lines = f.readlines()

names  = [x.split(',')[0] for x in lines[1:]]
hrrs   = [float(x.split(',')[1]) for x in lines[1:]]
angles = [float(x.split(',')[2].replace('\n','')) for x in lines[1:]]

hght_list = [x / 100 for x in range(5, 95, 10)]  #board height

hrrpuas = [x/(0.3*0.3) for x in hrrs]
#hrrpua = 555

templateFile = lines[0].split(',')[0]
with open(templateFile, 'r') as f:
    baseText = f.read()

for name, hrrpua, angle in zip(names, hrrpuas, angles):
    fileText = str(baseText)
    fileText = fileText.replace('ROTATION_ANGLE','ROTATION_ANGLE=%0.1f'%(angle))
    fileText = fileText.replace('HRRPUA','HRRPUA=%0.1f'%(hrrpua))
    fileText = fileText.replace('CHID','CHID=%s'%(name.replace('.fds','')))
    device_lines = generate_device_lines(angle, hght_list)
    fileText = fileText + '\n'.join(device_lines)
    
    with open("..//%s"%(name), 'w') as f:
        f.write(fileText)
