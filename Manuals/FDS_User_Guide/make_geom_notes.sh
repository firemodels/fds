#!/bin/bash

# Add LaTeX search path; Paths are ':' separated
export TEXINPUTS=".:../LaTeX_Style_Files:"

clean_build=1

# Build geom_notes

gitrevision=`git describe --long --dirty`
echo "\\newcommand{\\gitrevision}{$gitrevision}" > ../Bibliography/gitrevision.tex

pdflatex -interaction nonstopmode geom_notes &> geom_notes.err
bibtex geom_notes &> geom_notes.err
pdflatex -interaction nonstopmode geom_notes &> geom_notes.err
pdflatex -interaction nonstopmode geom_notes &> geom_notes.err

# make sure the guide exists
if [ ! -e geom_notes.pdf ]; then
  clean_build=0
  echo "***error: the geometry notes to build!"
fi

# Scan and report any errors in the LaTeX build process
if [[ `grep -E "Undefined control sequence|Error:|Fatal error|! LaTeX Error:|Paragraph ended before|Missing \\\$ inserted|Misplaced" -I geom_notes.err | grep -v "xpdf supports version 1.5"` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX errors detected:"
      grep -A 1 -E "Undefined control sequence|Error:|Fatal error|! LaTeX Error:|Paragraph ended before|Missing \\\$ inserted|Misplaced" -I geom_notes.err | grep -v "xpdf supports version 1.5"
      clean_build=0
fi

# Check for LaTeX warnings (undefined references or duplicate labels)
if [[ `grep -E "undefined|multiply defined|multiply-defined" -I geom_notes.err` == "" ]]
   then
      # Continue along
      :
   else
      echo "LaTeX warnings detected:"
      grep -E "undefined|multiply defined|multiply-defined" -I geom_notes.err
      clean_build=0
fi

if [[ $clean_build == 0 ]]
   then
      :
   else
      echo "geometry notes built successfully!"
fi    
