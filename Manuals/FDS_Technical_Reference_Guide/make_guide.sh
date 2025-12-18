#!/bin/bash

# Add LaTeX search path; Paths are ':' separated
export TEXINPUTS=".:../LaTeX_Style_Files:"

clean_build=1

# Build FDS Technical Reference Guide

gitrevision=`git describe --abbrev=7 --long --dirty`
echo "\\newcommand{\\gitrevision}{$gitrevision}" > ../Bibliography/gitrevision.tex

pdflatex -interaction nonstopmode FDS_Technical_Reference_Guide &> FDS_Technical_Reference_Guide.err
biber FDS_Technical_Reference_Guide &> FDS_Technical_Reference_Guide.err
pdflatex -interaction nonstopmode FDS_Technical_Reference_Guide &> FDS_Technical_Reference_Guide.err
pdflatex -interaction nonstopmode FDS_Technical_Reference_Guide &> FDS_Technical_Reference_Guide.err
pdflatex -interaction nonstopmode FDS_Technical_Reference_Guide &> FDS_Technical_Reference_Guide.err

# make sure the guide exists
if [ ! -e FDS_Technical_Reference_Guide.pdf ]; then
  clean_build=0
  echo "***error: the FDS Technical Reference Guide failed to build!"
fi

# Scan and report any errors in the LaTeX build process
if [[ `grep -E "Too many|Undefined control sequence|Error:|Fatal error|! LaTeX Error:|Paragraph ended before|Missing \\\$ inserted|Misplaced" -I FDS_Technical_Reference_Guide.err | grep -v "xpdf supports version 1.5"` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX errors detected:"
      grep -A 1 -E "Too many|Undefined control sequence|Error:|Fatal error|! LaTeX Error:|Paragraph ended before|Missing \\\$ inserted|Misplaced" -I FDS_Technical_Reference_Guide.err | grep -v "xpdf supports version 1.5"
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

   echo
   lacheck $texfile.tex

   # for other options see http://tex.loria.fr/outils/ChkTeX.pdf
   # http://tex.stackexchange.com/questions/46740/is-there-a-way-to-make-chktex-ignore-tikz-code
   # Newer ChkTeX: http://www.nongnu.org/chktex/
   # This does not look through files in the order they are called. Close enough.
   echo
   for ff in *.tex ; do
      echo $ff
      chktex -I0 -n1 -v4 $ff
   done
fi
