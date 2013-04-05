#!/bin/bash

# Add LaTeX search path; Paths are ':' separated
export TEXINPUTS=".:../LaTeX_Style_Files:"

clean_build=1

# Build FDS Technical Reference Guide
pdflatex -interaction nonstopmode FDS_Technical_Reference_Guide &> FDS_Technical_Reference_Guide.err
bibtex FDS_Technical_Reference_Guide &> FDS_Technical_Reference_Guide.err
pdflatex -interaction nonstopmode FDS_Technical_Reference_Guide &> FDS_Technical_Reference_Guide.err
pdflatex -interaction nonstopmode FDS_Technical_Reference_Guide &> FDS_Technical_Reference_Guide.err

# Scan and report any errors in the LaTeX build process
if [[ `grep -E "Error:|Fatal error|! LaTeX Error:" -I FDS_Technical_Reference_Guide.err | grep -v "xpdf supports version 1.5"` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX errors detected:"
      grep -E "Error:|Fatal error|! LaTeX Error:" -I FDS_Technical_Reference_Guide.err | grep -v "xpdf supports version 1.5"
      clean_build=0
fi

# Check for LaTeX warnings (undefined references or duplicate labels)
if [[ `grep -E "undefined|multiply defined|multiply-defined" -I FDS_Technical_Reference_Guide.err` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX warnings detected:"
      grep -E "new file|undefined|multiply defined|multiply-defined" -I FDS_Technical_Reference_Guide.err
      clean_build=0
fi

if [[ $clean_build == 0 ]]
   then
      :
   else
      echo "FDS Technical Reference Guide built successfully!"
fi    

if [ "$1" == "verbose" ]; then
   echo
   echo 'Verbose list of warnings:'
   echo
   texfile='FDS_Technical_Reference_Guide'
   grep -i error $texfile.log
   grep -i warning $texfile.log
   grep -i undefined $texfile.log
   grep 'not available' $texfile.log
   grep substituted $texfile.log
   grep 'not found' $texfile.log
   grep 'Rerun' $texfile.log
   grep 'Missing character' $texfile.log
   grep 'undefined' $texfile.log
   grep 'multiply defined' $texfile.log
   grep 'multiply-defined' $texfile.log
   grep 'Overfull \\hbox' $texfile.log
   grep 'Underfull \\hbox' $texfile.log
fi
