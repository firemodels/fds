#!/bin/bash

# Add LaTeX search path; Paths are ':' separated
export TEXINPUTS=".:../LaTeX_Style_Files:"
AUXUSER=../FDS_User_Guide/FDS_User_Guide.aux
PDFUSER=../FDS_User_Guide/FDS_User_Guide.pdf

if [ -e "$AUXUSER" ]; then
  cp $AUXUSER .
else
  echo "***warning: $AUXUSER does not exist. Build the FDS User's"
  echo "            Guide before building the Validation Guide"
fi

clean_build=1

# Build FDS Validation Guide

gitrevision=`git describe --abbrev=7 --long --dirty`
echo "\\newcommand{\\gitrevision}{$gitrevision}" > ../Bibliography/gitrevision.tex

pdflatex -interaction nonstopmode FDS_Validation_Guide &> FDS_Validation_Guide.err
biber                             FDS_Validation_Guide &> FDS_Validation_Guide_biber.err
pdflatex -interaction nonstopmode FDS_Validation_Guide &> FDS_Validation_Guide.err
pdflatex -interaction nonstopmode FDS_Validation_Guide &> FDS_Validation_Guide.err
pdflatex -interaction nonstopmode FDS_Validation_Guide &> FDS_Validation_Guide.err
cat FDS_Validation_Guide_biber.err >> FDS_Validation_Guide.err

# make sure the guide exists
if [ ! -e FDS_Validation_Guide.pdf ]; then
  clean_build=0
  echo "***error: the FDS Validation Guide failed to build!"
fi

if [ -e "$PDFUSER" ]; then
  cp $PDFUSER .
fi

# Scan and report any errors in the LaTeX build process
if [[ `grep -E "Too many|Undefined control sequence|Error:|Fatal error|! LaTeX Error:|Paragraph ended before|Missing \\\$ inserted|Misplaced" -I FDS_Validation_Guide.err | grep -v "xpdf supports version 1.5"` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX errors detected:"
      grep -A 1 -E "Too many|Undefined control sequence|Error:|Fatal error|! LaTeX Error:|Paragraph ended before|Missing \\\$ inserted|Misplaced" -I FDS_Validation_Guide.err | grep -v "xpdf supports version 1.5"
      clean_build=0
fi

# Check for LaTeX warnings (undefined references or duplicate labels)
if [[ `grep -E "undefined|WARNING|ERROR|multiply defined" -I FDS_Validation_Guide.err | grep -v RF1 | grep -v RF2 | grep -v LastPage` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX warnings detected:"
      grep -E "undefined|WARNING|ERROR|multiply defined" -I FDS_Validation_Guide.err | grep -v RF1 | grep -v RF2 | grep -v LastPage
      clean_build=0
fi

if [[ $clean_build == 0 ]]
   then
      :
   else
      echo "FDS Validation Guide built successfully!"
fi
