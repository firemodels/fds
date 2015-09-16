#!/bin/bash

# validation_git_stats.sh
# Kristopher Overholt
# 1/15/2014

# This script outputs a LaTeX file with a table of the FDS validation
# sets and their corresponding GIT information (i.e., when the FDS
# output files were last commited to the repository). This table
# is then included in the FDS Validation Guide.

CURRENT_DIR=`pwd`
if [ "$FDSSMV" == "" ] ; then
   FDSSMV=~/FDS-SMVgitclean
fi
cd $FDSSMV/Utilities/Scripts

# Name and location of output .tex file with validation GIT statistics
OUTPUT_TEX_FILE=$FDSSMV/Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/ScatterPlots/validation_git_stats.tex

# Table header
echo "\begin{longtable}[c]{|l|c|c|c|}" > $OUTPUT_TEX_FILE
echo "\caption[Validation GIT Statistics]{Validation GIT statistics for all data sets}" >> $OUTPUT_TEX_FILE
echo "\label{validation_git_stats}" >> $OUTPUT_TEX_FILE
echo "\\\\ \hline" >> $OUTPUT_TEX_FILE
echo "Dataset  &  FDS Revision Date  &  FDS Revision String&time (s)\\\\ \hline \hline" >> $OUTPUT_TEX_FILE
echo "\endfirsthead" >> $OUTPUT_TEX_FILE
echo "\hline" >> $OUTPUT_TEX_FILE
echo "Dataset  &  FDS Revision Date  &  FDS Revision String&time (s)\\\\ \hline \hline" >> $OUTPUT_TEX_FILE
echo "\endhead" >> $OUTPUT_TEX_FILE

./makegittable.sh >> $OUTPUT_TEX_FILE

# Table footer
echo "\end{longtable}" >> $OUTPUT_TEX_FILE

cd $CURRENT_DIR

