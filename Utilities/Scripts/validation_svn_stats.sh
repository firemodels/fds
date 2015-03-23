#!/bin/bash

# validation_svn_stats.sh
# Kristopher Overholt
# 1/15/2014

# This script outputs a LaTeX file with a table of the FDS validation
# sets and their corresponding SVN information (i.e., when the FDS
# output files were last commited to the repository). This table
# is then included in the FDS Validation Guide.

CURRENT_DIR=`pwd`
SVNROOT=~/FDS-SMV
cd $SVNROOT/Utilities/Scripts

# Name and location of output .tex file with validation SVN statistics
OUTPUT_TEX_FILE=$SVNROOT/Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/ScatterPlots/validation_svn_stats.tex

# Table header
echo "\begin{longtable}[c]{|l|c|c|}" > $OUTPUT_TEX_FILE
echo "\caption[Validation SVN Statistics]{Validation SVN statistics for all data sets}" >> $OUTPUT_TEX_FILE
echo "\label{validation_svn_stats}" >> $OUTPUT_TEX_FILE
echo "\\\\ \hline" >> $OUTPUT_TEX_FILE
echo "Dataset  &  FDS Build Date  &  FDS Revision \\\\ \hline \hline" >> $OUTPUT_TEX_FILE
echo "\endfirsthead" >> $OUTPUT_TEX_FILE
echo "\hline" >> $OUTPUT_TEX_FILE
echo "Dataset  &  FDS Build Date  &  FDS Revision \\\\ \hline \hline" >> $OUTPUT_TEX_FILE
echo "\endhead" >> $OUTPUT_TEX_FILE

./makesvntable.sh >> $OUTPUT_TEX_FILE

# Table footer
echo "\end{longtable}" >> $OUTPUT_TEX_FILE

cd $CURRENT_DIR

