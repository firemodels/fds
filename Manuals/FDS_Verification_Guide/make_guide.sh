#!/bin/bash

# Add LaTeX search path; Paths are ':' separated
export TEXINPUTS=".:../LaTeX_Style_Files:"
AUXUSER=../FDS_User_Guide/FDS_User_Guide.aux
PDFVER=FDS_Verification_Guide.pdf

if [ -e "$AUXUSER" ]; then
  cp $AUXUSER .
else
  echo "***warning: $AUXUSER does not exist. Build the FDS User's"
  echo "            guide before building the Verification guide"
fi

clean_build=1

# Build FDS Verification Guide

gitrevision=`git describe --long --dirty`
echo "\\newcommand{\\gitrevision}{$gitrevision}" > ../Bibliography/gitrevision.tex

pdflatex -interaction nonstopmode FDS_Verification_Guide &> FDS_Verification_Guide.err
bibtex FDS_Verification_Guide &> FDS_Verification_Guide.err
pdflatex -interaction nonstopmode FDS_Verification_Guide &> FDS_Verification_Guide.err
pdflatex -interaction nonstopmode FDS_Verification_Guide &> FDS_Verification_Guide.err
if [ -e "$PDFVER" ]; then
  cp $PDFVER ../FDS_User_Guide/.
fi

# Scan and report any errors in the LaTeX build process
if [[ `grep -E "Error:|Fatal error|! LaTeX Error:|Paragraph ended before|Missing \\\$ inserted|Misplaced" -I FDS_Verification_Guide.err | grep -v "xpdf supports version 1.5"` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX errors detected:"
      grep -E "Error:|Fatal error|! LaTeX Error:|Paragraph ended before|Missing \\\$ inserted|Misplaced" -I FDS_Verification_Guide.err | grep -v "xpdf supports version 1.5"
      clean_build=0
fi

# Check for LaTeX warnings (undefined references or duplicate labels)
if [[ `grep -E "undefined|multiply defined" -I FDS_Verification_Guide.err | grep -v RF1 | grep -v LastPage` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX warnings detected:"
      grep -E "undefined|multiply defined" -I FDS_Verification_Guide.err | grep -v RF1 | grep -v LastPage
      clean_build=0
fi

if [[ $clean_build == 0 ]]
   then
      :
   else
      echo "FDS Verification Guide built successfully!"
fi    
