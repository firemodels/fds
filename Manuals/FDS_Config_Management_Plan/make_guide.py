
import os
import subprocess
import re

os.environ["TEXINPUTS"] = ".:../LaTeX_Style_Files:"

clean_build = 1

try:
    gitrevision = subprocess.check_output(
        ["git", "describe", "--abbrev=7", "--long", "--dirty"],
        stderr=subprocess.STDOUT,
        text=True
    ).strip()
except subprocess.CalledProcessError:
    gitrevision = ""

with open("../Bibliography/gitrevision.tex", "w") as f:
    f.write(f"\\newcommand{{\\gitrevision}}{{{gitrevision}}}\n")

with open("FDS_Config_Management_Plan.err", "w") as f_err:
    subprocess.run(
        ["pdflatex", "-interaction", "nonstopmode", "FDS_Config_Management_Plan"],
        stdout=f_err,
        stderr=f_err
    )

with open("FDS_Config_Management_Plan_biber.err", "w") as f_biber_err:
    subprocess.run(
        ["biber", "FDS_Config_Management_Plan"],
        stdout=f_biber_err,
        stderr=f_biber_err
    )

with open("FDS_Config_Management_Plan.err", "w") as f_err:
    subprocess.run(
        ["pdflatex", "-interaction", "nonstopmode", "FDS_Config_Management_Plan"],
        stdout=f_err,
        stderr=f_err
    )

with open("FDS_Config_Management_Plan.err", "w") as f_err:
    subprocess.run(
        ["pdflatex", "-interaction", "nonstopmode", "FDS_Config_Management_Plan"],
        stdout=f_err,
        stderr=f_err
    )

with open("FDS_Config_Management_Plan_biber.err", "r") as f_src:
    content = f_src.read()
with open("FDS_Config_Management_Plan.err", "a") as f_dst:
    f_dst.write(content)

if not os.path.exists("FDS_Config_Management_Plan.pdf"):
    clean_build = 0
    print("***error: the FDS Config Management Plan failed to build!")

error_regex = re.compile(r"Too many|Undefined control sequence|Error:|Fatal error|! LaTeX Error:|Paragraph ended before|Missing \$ inserted|Misplaced")
exclusion = "xpdf supports version 1.5"

def check_for_errors(filename):
    try:
        with open(filename, "r", errors="ignore") as f:
            for line in f:
                if error_regex.search(line):
                    if exclusion not in line:
                        return True
    except FileNotFoundError:
        pass
    return False

if not check_for_errors("FDS_Config_Management_Plan.err"):
    # Continue along
    pass
else:
    print("LaTeX errors detected:")
    try:
        with open("FDS_Config_Management_Plan.err", "r", errors="ignore") as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if error_regex.search(line):
                    context_lines = lines[i : i + 2]
                    for cl in context_lines:
                        if exclusion not in cl:
                            print(cl.rstrip())
    except FileNotFoundError:
        pass
    clean_build = 0

warning_regex = re.compile(r"undefined|WARNING|ERROR|multiply defined|multiply-defined")

def check_for_warnings(filename):
    matches = []
    try:
        with open(filename, "r", errors="ignore") as f:
            for line in f:
                if warning_regex.search(line):
                    matches.append(line.rstrip())
    except FileNotFoundError:
        pass
    return matches

warnings_found = check_for_warnings("FDS_Config_Management_Plan.err")
if not warnings_found:
    # Continue along
    pass
else:
    print("LaTeX warnings detected:")
    for warn in warnings_found:
        print(warn)
    clean_build = 0

if clean_build == 0:
    pass
else:
    print("FDS Config Management Plan built successfully!")

