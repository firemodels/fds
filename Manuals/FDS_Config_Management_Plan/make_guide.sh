#!/bin/bash

# Add LaTeX search path; Paths are ':' separated
export TEXINPUTS=".:../LaTeX_Style_Files:"

clean_build=1

# Build FDS Config Management Plan

gitrevision=`git describe --abbrev=7 --long --dirty`
echo "\\newcommand{\\gitrevision}{$gitrevision}" > ../Bibliography/gitrevision.tex

pdflatex -interaction nonstopmode FDS_Config_Management_Plan &> FDS_Config_Management_Plan.err
biber                             FDS_Config_Management_Plan &> FDS_Config_Management_Plan_biber.err
pdflatex -interaction nonstopmode FDS_Config_Management_Plan &> FDS_Config_Management_Plan.err
pdflatex -interaction nonstopmode FDS_Config_Management_Plan &> FDS_Config_Management_Plan.err
pdflatex -interaction nonstopmode FDS_Config_Management_Plan &> FDS_Config_Management_Plan.err
cat FDS_Config_Management_Plan_biber.err >> FDS_Config_Management_Plan.err

# make sure the guide exists
if [ ! -e FDS_Config_Management_Plan.pdf ]; then
  clean_build=0
  echo "***error: the FDS Config Management Plan failed to build!"
fi
  
# Scan and report any errors in the LaTeX build process
if [[ `grep -E "Too many|Undefined control sequence|Error:|Fatal error|! LaTeX Error:|Paragraph ended before|Missing \\\$ inserted|Misplaced" -I FDS_Config_Management_Plan.err | grep -v "xpdf supports version 1.5"` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX errors detected:"
      grep -A 1 -E "Too many|Undefined control sequence|Error:|Fatal error|! LaTeX Error:|Paragraph ended before|Missing \\\$ inserted|Misplaced" -I FDS_Config_Management_Plan.err | grep -v "xpdf supports version 1.5"
      clean_build=0
fi

# Check for LaTeX warnings (undefined references or duplicate labels)
if [[ `grep -E "undefined|WARNING|ERROR|multiply defined|multiply-defined" -I FDS_Config_Management_Plan.err` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX warnings detected:"
      grep -E "undefined|WARNING|ERROR|multiply defined|multiply-defined" -I FDS_Config_Management_Plan.err
      clean_build=0
fi

if [[ $clean_build == 0 ]]
   then
      :
   else
      echo "FDS Config Management Plan built successfully!"
fi
