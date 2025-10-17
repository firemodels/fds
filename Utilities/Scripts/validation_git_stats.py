#!/usr/bin/env python3

import os
import sys
import subprocess
import re
import tempfile
from pathlib import Path
import argparse

def MAKEGITENTRY(DIR, FIREMODELS_ROOT, TEMPDIR):
    """
    This function outputs a LaTeX table entry with git information for a validation set.
    """
    gitrevisions = os.path.join(TEMPDIR, f'gitrevisions.{os.getpid()}')

    # Collect all git.txt files and sort uniquely
    git_txt_pattern = os.path.join(FIREMODELS_ROOT, 'out', DIR, '*git.txt')

    # Use shell to expand glob and collect content
    try:
        cmd = f'cat {git_txt_pattern} 2> /dev/null | sort -u'
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        with open(gitrevisions, 'w') as f:
            f.write(result.stdout)
    except Exception:
        # If no files found, create empty file
        with open(gitrevisions, 'w') as f:
            f.write('')

    # Read first line
    gitrev = ''
    try:
        with open(gitrevisions, 'r') as f:
            gitrev = f.readline().strip()
    except Exception:
        gitrev = ''

    output = ''
    if gitrev != '':
        # Extract git revision short hash
        # awk -F - '{print $(NF-1)}' | sed 's/^g//'
        parts = gitrev.split('-')
        if len(parts) >= 2:
            gitrevshort = parts[-2]
            gitrevshort = gitrevshort[1:] if gitrevshort.startswith('g') else gitrevshort
        else:
            gitrevshort = gitrev

        # Get git date
        gitdate = ''
        gitdate2 = ''
        try:
            result = subprocess.run(
                ['git', 'show', '-s', '--format=%aD', gitrevshort],
                capture_output=True,
                text=True,
                stderr=subprocess.DEVNULL
            )
            if result.returncode == 0 and result.stdout.strip():
                # Parse date: awk '{print $3,$2",",$4}'
                date_parts = result.stdout.strip().split()
                if len(date_parts) >= 5:
                    gitdate = f"{date_parts[2]} {date_parts[1]}, {date_parts[3]}"
        except Exception:
            gitdate = ''

        if gitdate == '':
            gitdate = 'undefined'
            gitdate2 = '2000000000'

            # Check if ~/FDS-SMV exists
            fds_smv_path = os.path.expanduser('~/FDS-SMV')
            if os.path.exists(fds_smv_path):
                CUR_DIR = os.getcwd()
                os.chdir(fds_smv_path)

                # Extract different part of gitrev: awk -F - '{print $4}' | sed 's/^.\{1\}//'
                parts = gitrev.split('-')
                if len(parts) >= 4:
                    gitrevshort = parts[3]
                    gitrevshort = gitrevshort[1:] if len(gitrevshort) > 0 else gitrevshort

                gitdateold = ''
                try:
                    result = subprocess.run(
                        ['git', 'show', '-s', '--format=%aD', gitrevshort],
                        capture_output=True,
                        text=True,
                        stderr=subprocess.DEVNULL
                    )
                    if result.returncode == 0 and result.stdout.strip():
                        date_parts = result.stdout.strip().split()
                        if len(date_parts) >= 5:
                            gitdateold = f"{date_parts[2]} {date_parts[1]}, {date_parts[3]}"
                except Exception:
                    gitdateold = ''

                if gitdateold != '':
                    gitdate = gitdateold
                    try:
                        result = subprocess.run(
                            ['git', 'show', '-s', '--format=%at', gitrevshort],
                            capture_output=True,
                            text=True
                        )
                        if result.returncode == 0 and result.stdout.strip():
                            gitdate2 = result.stdout.strip().split()[0]
                    except Exception:
                        pass

                os.chdir(CUR_DIR)
        else:
            # Get Unix timestamp
            try:
                result = subprocess.run(
                    ['git', 'show', '-s', '--format=%at', gitrevshort],
                    capture_output=True,
                    text=True
                )
                if result.returncode == 0 and result.stdout.strip():
                    gitdate2 = result.stdout.strip().split()[0]
            except Exception:
                gitdate2 = ''

        # Escape underscores for LaTeX
        dir_escaped = DIR.replace('_', '\\_')
        output = f"{dir_escaped}  & {gitdate} & {gitrev} & {gitdate2} \\\\ \\hline\n"

    # Remove temporary file
    try:
        os.remove(gitrevisions)
    except Exception:
        pass

    return output


def main():
    """
    Main function that generates the LaTeX table with validation git statistics.
    """
    CURRENT_DIR = os.getcwd()

    # Determine repo root
    SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))
    os.chdir(os.path.join(SCRIPTDIR, '../../..'))
    FIREMODELS_ROOT = os.getcwd()

    # Set up temp directory
    TEMPDIR = os.path.join(os.path.expanduser('~'), 'temp')
    if not os.path.exists(TEMPDIR):
        os.makedirs(TEMPDIR)

    # Set environment variable
    os.environ['FIREMODELS_ROOT'] = FIREMODELS_ROOT

    # Parse command line options (kept for backwards compatibility)
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', dest='IGNORE', help='Ignored parameter for backwards compatibility')
    args = parser.parse_args()

    # Change to scripts directory
    os.chdir(os.path.join(FIREMODELS_ROOT, 'fds/Utilities/Scripts'))

    # Name and location of output .tex file with validation GIT statistics
    OUTPUT_TEX_FILE = os.path.join(
        FIREMODELS_ROOT,
        'fds/Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/ScatterPlots/validation_git_stats.tex'
    )

    # Table header
    with open(OUTPUT_TEX_FILE, 'w') as f:
        f.write("\\begin{longtable}[c]{|l|c|c|}\n")
        f.write("\\caption[Validation Git Statistics]{Validation Git statistics for all data sets}\n")
        f.write("\\label{validation_git_stats}\n")
        f.write("\\\\ \\hline\n")
        f.write("Dataset  &  FDS Revision Date  &  FDS Revision String\\\\ \\hline \\hline\n")
        f.write("\\endfirsthead\n")
        f.write("\\hline\n")
        f.write("Dataset  &  FDS Revision Date  &  FDS Revision String\\\\ \\hline \\hline\n")
        f.write("\\endhead\n")

    # Table body
    maketable = os.path.join(FIREMODELS_ROOT, 'fds/Validation/Process_All_Output.sh')
    CASELIST = os.path.join(TEMPDIR, f'temp.out.{os.getpid()}')
    TABLE_ENTRIES = os.path.join(TEMPDIR, f'temp2.out.{os.getpid()}')

    # Extract case list from Process_All_Output.sh
    # grep PROCESS $maketable | awk 'BEGIN { FS = " " } ; { print $2 }' | awk '{if(NR>1)print}'
    try:
        with open(maketable, 'r') as f:
            lines = f.readlines()

        cases = []
        line_num = 0
        for line in lines:
            if 'PROCESS' in line:
                parts = line.strip().split()
                if len(parts) >= 2:
                    line_num += 1
                    if line_num > 1:  # Skip first match
                        cases.append(parts[1])

        with open(CASELIST, 'w') as f:
            for case in cases:
                f.write(f"{case}\n")
    except Exception as e:
        # If file doesn't exist or error, create empty caselist
        with open(CASELIST, 'w') as f:
            f.write('')

    # Process each case and generate table entries
    with open(TABLE_ENTRIES, 'w') as outf:
        try:
            with open(CASELIST, 'r') as inf:
                for line in inf:
                    p = line.strip()
                    if p:
                        entry = MAKEGITENTRY(p, FIREMODELS_ROOT, TEMPDIR)
                        if entry:
                            outf.write(entry)
        except Exception:
            pass

    # Sort table entries and append to output file
    # cat $TABLE_ENTRIES | sort -n -t '&' -k 4 | awk -F "&" '{ print $1 "&" $2 "&" $3 "\\\\ \\hline"}'
    try:
        with open(TABLE_ENTRIES, 'r') as f:
            entries = f.readlines()

        # Sort by 4th field (numeric, using & as delimiter)
        def sort_key(line):
            parts = line.split('&')
            if len(parts) >= 4:
                # Extract numeric value from 4th field
                try:
                    return int(parts[3].strip().split()[0])
                except:
                    return 0
            return 0

        sorted_entries = sorted(entries, key=sort_key)

        with open(OUTPUT_TEX_FILE, 'a') as f:
            for entry in sorted_entries:
                parts = entry.split('&')
                if len(parts) >= 4:
                    # Reconstruct line with first 3 fields
                    output_line = f"{parts[0]}&{parts[1]}&{parts[2]}\\\\ \\hline\n"
                    f.write(output_line)
    except Exception:
        pass

    # Clean up temporary files
    try:
        os.remove(CASELIST)
    except Exception:
        pass
    try:
        os.remove(TABLE_ENTRIES)
    except Exception:
        pass

    # Table footer
    with open(OUTPUT_TEX_FILE, 'a') as f:
        f.write("\\end{longtable}\n")

    # Return to original directory
    os.chdir(CURRENT_DIR)


if __name__ == '__main__':
    main()

