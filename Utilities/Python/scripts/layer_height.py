import os
import numpy as np

validation_dir = '../../Validation/'
output_base_dir = '../../../out/'

list_dir = os.listdir(validation_dir)
output_dir = []
input_file = []
k = 0

for name in list_dir:
    output_directory = os.path.join(output_base_dir, name)
    if os.path.exists(output_directory) and os.path.isdir(output_directory):
        list_files = [f for f in os.listdir(output_directory) if f.endswith('HGL.input')]
        if list_files:
            for file_name in list_files:
                if not file_name.startswith('.'):  # ignore hidden files
                    k += 1
                    output_dir.append(os.path.join(output_base_dir, name) + '/')
                    input_file.append(file_name)

# Uncomment the following line to just list the files for testing purposes
# exit()

ntd = 20000
ncd = 500

for idx in range(len(input_file)):  # input_file loop
    tmp = np.zeros((ntd, ncd))
    
    input_path = os.path.join(output_dir[idx], input_file[idx])
    with open(input_path, 'r') as fid:
        # Read number of trees
        ntrees = int(fid.readline())
        
        # Read number of TCs in the tree
        ntc = int(fid.readline())
        ztc = np.zeros(ntc)
        icol = np.zeros((ntc, ntrees), dtype=int)
        for n in range(ntc):
            line = fid.readline()
            S = line.strip().split()
            ztc[n] = float(S[0])
            for nn in range(ntrees):
                icol[n, nn] = int(S[nn + 1]) - 1
        
        # Read weight of each tree
        line = fid.readline()
        S = line.strip().split()
        wgt = np.array([float(s) for s in S[:ntrees]])
        
        # Read data file name
        infile = fid.readline().strip()
        
        # Read number of columns in data file
        nc = int(fid.readline())
        
        # Read row number where data starts
        nr = int(fid.readline())
        
        # Read ceiling height
        z_h = float(fid.readline())
        
        # Read starting time
        t_start = float(fid.readline())
        
        # Read name of output file
        outfile = fid.readline().strip()
    
    # Read data from file
    data_path = os.path.join(output_dir[idx], infile)
    M = np.loadtxt(data_path, delimiter=',', skiprows=nr-1)
    t = M[:, 0]
    d = M[:, 1:nc]
    
    z = np.zeros(ntc)
    for n in range(ntc - 1):
        z[n] = (ztc[n] + ztc[n + 1]) / 2
    z[-1] = z_h
    
    fout_path = os.path.join(output_dir[idx], outfile)
    with open(fout_path, 'w') as fout:
        fout.write('Time, Height, T_lower, T_upper\n')
        
        for i in range(len(t)):  # time loop
            if t[i] < t_start:
                continue
            tmp[i, :] = 0
            for nn in range(ntrees):
                for n in range(ntc):
                    tmp[i, n] += (273 + d[i, icol[n, nn] - 1]) * wgt[nn]
            
            i1 = 0
            i2 = 0
            for n in range(ntc):
                if n == 0:
                    i1 += tmp[i, n] * (z[n] - 0)
                    i2 += (1.0 / tmp[i, n]) * (z[n] - 0)
                else:
                    i1 += tmp[i, n] * (z[n] - z[n - 1])
                    i2 += (1.0 / tmp[i, n]) * (z[n] - z[n - 1])
            zint = tmp[i, 0] * (i1 * i2 - z_h**2) / (0.00001 + i1 + i2 * tmp[i, 0]**2 - 2 * tmp[i, 0] * z_h)
            tmpl = tmp[i, 0]
            i1 = 0
            for n in range(ntc):
                if n == 0:
                    if z[n] > zint:
                        if 0 >= zint:
                            i1 += tmp[i, n] * (z[n] - 0)
                        if 0 < zint:
                            i1 += tmp[i, n] * (z[n] - zint)
                else:
                    if z[n] > zint:
                        if z[n - 1] >= zint:
                            i1 += tmp[i, n] * (z[n] - z[n - 1])
                        if z[n - 1] < zint:
                            i1 += tmp[i, n] * (z[n] - zint)
            tmph = max(tmpl, (1.0 / (z_h - zint)) * i1)
            
            fout.write(f'{t[i]}, {zint}, {tmpl - 273}, {tmph - 273}\n')

