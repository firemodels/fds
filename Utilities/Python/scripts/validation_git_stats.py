
# Generate the LaTeX table with validation git statistics.

import subprocess
from pathlib import Path
import glob

outdir = '../../../out/'
valdir = '../../Validation/'
resdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/ScatterPlots/'

def MAKEGITENTRY(case_name):

    # Output a single LaTeX table entry with git information for a validation set.

    # Collect all git.txt files and sort uniquely
    git_file_pattern = outdir + case_name + '/*git.txt'
    matching_files = set(Path().glob(git_file_pattern))

    # Read first line
    gitrev = ''
    for file_path in sorted(matching_files):
        with open(file_path.as_posix(), 'r') as fff:
            gitrev = fff.readline().strip()

    output = ''
    if gitrev != '':
        # Extract git revision short hash
        parts = gitrev.split('-')
        if len(parts) >= 2:
            gitrevshort = parts[-2]
            gitrevshort = gitrevshort[1:] if gitrevshort.startswith('g') else gitrevshort
        else:
            gitrevshort = gitrev

        # Get git date
        gitdate = ''
        try:
            result = subprocess.run( ['git', 'show', '-s', '--format=%aD', gitrevshort], capture_output=True, text=True)
            if result.returncode == 0 and result.stdout.strip():
                date_parts = result.stdout.strip().split()
                if len(date_parts) >= 5:
                    gitdate = f"{date_parts[2]} {date_parts[1]}, {date_parts[3]}"
        except Exception:
            gitdate = 'Unknown'

        # Escape underscores for LaTeX
        dir_escaped = case_name.replace('_', '\\_')
        output = f"{dir_escaped}  & {gitdate} & {gitrev} \\\\ \\hline\n"

    return output


# Create a LaTeX table

OUTPUT_TEX_FILE = resdir + 'validation_git_stats.tex'

with open(OUTPUT_TEX_FILE, 'w') as outf:
    outf.write("\\begin{longtable}[c]{|l|c|c|}\n")
    outf.write("\\caption[Validation Git Statistics]{Validation Git statistics for all data sets}\n")
    outf.write("\\label{validation_git_stats}\n")
    outf.write("\\\\ \\hline\n")
    outf.write("Dataset  &  FDS Revision Date  &  FDS Revision String\\\\ \\hline \\hline\n")
    outf.write("\\endfirsthead\n")
    outf.write("\\hline\n")
    outf.write("Dataset  &  FDS Revision Date  &  FDS Revision String\\\\ \\hline \\hline\n")
    outf.write("\\endhead\n")


# Extract case list from Validation/Process_All_Output.sh

with open(valdir + 'Process_All_Output.sh', 'r') as inf:
    lines = inf.readlines()

cases = []
line_num = 0
for line in lines:
    if 'PROCESS' in line:
        parts = line.strip().split()
        if len(parts) >= 2:
            line_num += 1
            if line_num > 1:  # Skip first match
                cases.append(parts[1])

# Process each case and generate table entries

for case in cases:
    entry = MAKEGITENTRY(case)
    with open(OUTPUT_TEX_FILE, 'a') as outf:
        outf.write(entry)

# Table footer

with open(OUTPUT_TEX_FILE, 'a') as f:
    f.write("\\end{longtable}\n")

