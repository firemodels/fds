
# Generate the LaTeX table with validation git statistics and compare scatterplot
# statistics with the version-controlled baseline.

import csv
import subprocess
import sys
from pathlib import Path
import glob
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation

outdir = '../../../out/'
valdir = '../../Validation/'
resdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Scatterplots/'


def compare_scatterplot_statistics():
    """Report changes over 10% as errors only when model agreement worsens."""
    output_file = Path(resdir) / 'validation_scatterplot_output.csv'
    baseline_file = Path(resdir) / 'validation_scatterplot_output_baseline.csv'
    metrics = ('Sigma_Model', 'Bias')
    difference_columns = [metric + '_Relative_Difference' for metric in metrics]

    try:
        with baseline_file.open(newline='') as inf:
            baseline = {row['Quantity']: row for row in csv.DictReader(inf)}
        with output_file.open(newline='') as inf:
            reader = csv.DictReader(inf)
            fieldnames = list(reader.fieldnames or [])
            rows = list(reader)
    except OSError as exc:
        print(f'Error: validation_git_stats: cannot compare scatterplot statistics: {exc}')
        return

    # Re-running the script updates the existing columns instead of duplicating them.
    fieldnames += [name for name in difference_columns if name not in fieldnames]
    for row in rows:
        quantity = row['Quantity']
        for column in difference_columns:
            row[column] = ''
        if quantity not in baseline:
            print(f'Error: validation_git_stats: {quantity}: missing scatterplot baseline.')
            continue

        for metric, column in zip(metrics, difference_columns):
            try:
                # Decimal keeps changes of exactly 10% from failing due to roundoff.
                current = Decimal(row[metric])
                reference = Decimal(baseline[quantity][metric])
                if not current.is_finite() or not reference.is_finite():
                    raise ValueError('non-finite statistic')
            except (KeyError, TypeError, ValueError, InvalidOperation):
                print(f'Error: validation_git_stats: {quantity}: invalid {metric} '
                      'in scatterplot output or baseline.')
                continue

            if reference == 0:
                difference = Decimal(0) if current == 0 else Decimal('Infinity').copy_sign(current)
            else:
                difference = (current - reference) / abs(reference)
            row[column] = f'{difference:.6f}'
            if abs(difference) > Decimal('0.10'):
                # Lower model scatter is better; bias is better closer to one.
                # A bias crossing one can change significantly without worsening.
                current_error = current if metric == 'Sigma_Model' else abs(current - 1)
                baseline_error = reference if metric == 'Sigma_Model' else abs(reference - 1)
                severity = 'Error' if current_error > baseline_error else 'Warning'
                print(f'{severity}: validation_git_stats: {quantity}: {metric} relative '
                      f'difference {difference:.2%} exceeds 10% in magnitude '
                      f'(current={current}, baseline={reference}).')

    with output_file.open('w', newline='') as outf:
        writer = csv.DictWriter(outf, fieldnames=fieldnames, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


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
    gitdate = 'Unknown'
    # Keep unknown dates comparable with timezone-aware Git dates.
    git_datetime = datetime.min.replace(tzinfo=timezone.utc)

    if gitrev != '':
        # Extract git revision short hash
        parts = gitrev.split('-')
        if len(parts) >= 2:
            gitrevshort = parts[-2]
            gitrevshort = gitrevshort[1:] if gitrevshort.startswith('g') else gitrevshort
        else:
            gitrevshort = gitrev

        # Get git date
        try:
            result = subprocess.run(
                ['git', 'show', '-s', '--format=%aD', gitrevshort],
                capture_output=True,
                text=True
            )

            if result.returncode == 0 and result.stdout.strip():
                date_string = result.stdout.strip()

                # Convert Git's RFC 2822 date to a datetime object.
                # This is used for sorting, rather than sorting the
                # formatted date string alphabetically.
                git_datetime = datetime.strptime(
                    date_string,
                    '%a, %d %b %Y %H:%M:%S %z'
                )

                # Format date for LaTeX table
                gitdate = git_datetime.strftime('%B %d, %Y')

            else:
                reason = result.stderr.strip() or 'git show returned no date'
                print(f'[validation_git_stats] {case_name}: cannot determine date '
                      f'for revision {gitrevshort}: {reason}', file=sys.stderr)

        except Exception as exc:
            gitdate = 'Unknown'
            git_datetime = datetime.min.replace(tzinfo=timezone.utc)
            print(f'[validation_git_stats] {case_name}: cannot determine date '
                  f'for revision {gitrevshort}: {exc}', file=sys.stderr)

        # Escape underscores for LaTeX
        dir_escaped = case_name.replace('_', '\\_')

        output = f"{dir_escaped}  & {gitdate} & {gitrev} \\\\ \\hline\n"

    return git_datetime, output


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


# Process each case and collect entries

entries = []

for case in cases:
    git_datetime, entry = MAKEGITENTRY(case)

    if entry:
        entries.append((git_datetime, entry))


# Sort cases from oldest to newest by Git revision date

entries.sort(key=lambda x: x[0])


# Write sorted entries to the LaTeX table

with open(OUTPUT_TEX_FILE, 'a') as outf:
    for git_datetime, entry in entries:
        outf.write(entry)


# Table footer

with open(OUTPUT_TEX_FILE, 'a') as f:
    f.write("\\end{longtable}\n")


compare_scatterplot_statistics()
