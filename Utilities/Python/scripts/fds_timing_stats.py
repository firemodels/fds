
# Produce a table of verification case CPU times in Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/Scatterplots

import os
import subprocess
import csv

curdir = os.getcwd()
verdir = '../../Verification/'
resdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/Scatterplots/'
resfile = resdir + 'fds_timing_stats.csv'

with open(resfile, 'w') as f:
    f.write('FDS Case,Wall Clock Time (s),CPU Time (s),Number of Cells,Number of Time Steps,Performance Metric (1e-6)\n')

with open(verdir + 'FDS_Cases.sh', 'r') as casefile:
    for line in casefile:
        if line.strip().startswith('$QFDS'):
            parts = line.split()
            if '-d' in parts:
                d_index = parts.index('-d')
                if d_index + 2 < len(parts):
                    subdir = parts[d_index + 1]
                    fdsfile = parts[d_index + 2]

                    # Extract base filename without extension
                    fdsbase = verdir + subdir + '/' + os.path.splitext(fdsfile)[0]
                    outfile = f"{fdsbase}.out"
                    cpufile = f"{fdsbase}_cpu.csv"

                    # Check if outfile exists, try alternate name if not
                    if not os.path.exists(outfile):
                        outfile = f"{fdsbase}_cat.out"
                        if not os.path.exists(outfile):
                            continue

                    # Check if cpufile exists, try alternate name if not
                    if not os.path.exists(cpufile):
                        cpufile = f"{fdsbase}_cat_cpu.csv"
                        if not os.path.exists(cpufile):
                            continue
                    
                    # Grep for wall clock time
                    WALL_CLOCK_TIME_VALUE = ""
                    with open(outfile, 'r') as f:
                        for line in f:
                            if "Total Elapsed Wall Clock Time (s):" in line:
                                WALL_CLOCK_TIME_VALUE = line.strip().split()[-1]
                    
                    # Grep for CPU time and units
                    TOTAL_CPU_TIME = 0.0
                    CPU_TIME_VALUES = []
                    
                    with open(cpufile, 'r') as f:
                        lines = f.readlines()
                        # Skip first line (header), process rest
                        for line in lines[1:]:
                            if line.strip():
                                # Get last column (split by comma)
                                parts = line.strip().split(',')
                                if parts:
                                    CPU_TIME_VALUES.append(parts[-1])
                    
                    # Process each CPU time value
                    for j in CPU_TIME_VALUES:
                        TOTAL_CPU_TIME += eval(j)
                    
                    CPU_TIME = TOTAL_CPU_TIME
                    
                    # Grep for number of cells in each dimension
                    X_CELLS = []
                    Y_CELLS = []
                    Z_CELLS = []
                    
                    with open(outfile, 'r') as f:
                        for line in f:
                            if "Cells in the X" in line:
                                X_CELLS.append(line.strip().split()[-1])
                            elif "Cells in the Y" in line:
                                Y_CELLS.append(line.strip().split()[-1])
                            elif "Cells in the Z" in line:
                                Z_CELLS.append(line.strip().split()[-1])
                    
                    numx = len(X_CELLS)
                    
                    # Sum over the number of cells (for multi-mesh cases)
                    NUM_TOTAL_CELLS = 0
                    for i in range(numx):
                        XI = int(X_CELLS[i])
                        YI = int(Y_CELLS[i])
                        ZI = int(Z_CELLS[i])
                        sumxyz = XI * YI * ZI
                        NUM_TOTAL_CELLS += sumxyz
                    
                    # Grep for number of time steps
                    NUM_TIME_STEPS = 0
                    time_step_lines = []
                    with open(outfile, 'r') as f:
                        for line in f:
                            if "Time Step     " in line:
                                time_step_lines.append(line)
                    
                    if time_step_lines:
                        # Get the last occurrence
                        last_line = time_step_lines[-1]
                        parts = last_line.strip().split()
                        if len(parts) >= 5:
                            # Get the 5th element from the end (NF-4 in awk)
                            NUM_TIME_STEPS = int(parts[-5])
                    
                    # Calculate nondimensional performance metric
                    # Skip over cases with no time steps
                    if NUM_TIME_STEPS == 0:
                        NUM_TIME_STEPS = 0
                        PERFORMANCE = 0
                    else:
                        # Calculate performance metric
                        PERFORMANCE = int(1000000 * TOTAL_CPU_TIME / (NUM_TOTAL_CELLS * NUM_TIME_STEPS))
                    
                    # Write results to fds_timing_stats.csv file
                    with open(resfile, 'a', newline='') as f:
                        csv_writer = csv.writer(f)
                        csv_writer.writerow([fdsfile,WALL_CLOCK_TIME_VALUE,CPU_TIME,NUM_TOTAL_CELLS,NUM_TIME_STEPS,PERFORMANCE])


# Sum up the wall clock times in the second column

TOTAL_CPU_TIME = 0.0
with open(resfile, 'r') as f:
    lines = f.readlines()
    for line in lines[1:]:
        if line.strip():
            fields = line.split(',')
            if len(fields) >= 3:
                j = fields[1].strip()
                TOTAL_CPU_TIME = TOTAL_CPU_TIME + eval(j)

# Get git short commit hash and append to tmpout

git_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=curdir, universal_newlines=True).strip()

# Append git hash and total CPU time to fds_timing_stats.csv

with open(resfile, 'a') as f:
    f.write(git_hash + '\n')
    f.write(str(TOTAL_CPU_TIME) + '\n')

