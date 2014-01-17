#!/bin/bash

# validation_svn_stats.sh
# Kristopher Overholt
# 1/15/2014

# This script outputs a LaTeX file with a table of the FDS validation
# sets and their corresponding SVN information (i.e., when the FDS
# output files were last commited to the repository). This table
# is then included in the FDS Validation Guide.

CURRENT_DIR=`pwd`
SVNROOT=`pwd`/../..
cd $SVNROOT/Validation

# Name and location of output .tex file with validation SVN statistics
OUTPUT_TEX_FILE=$SVNROOT/Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/ScatterPlots/validation_svn_stats.tex

# Generate arrays with the name of validation data sets, SVN revision number, and SVN revision date
VALIDATION_SETS=(`grep '$VDIR' Process_All_Output.sh | grep -v "#" | xargs -n 1 dirname | xargs -n 1 dirname | xargs -n 1 basename | xargs -i svn info {}/FDS_Output_Files | awk '{if($0 != ""){ if(s){s=s"*"$0}else{s=$0}}else{ print s"*";s=""}}END{print s"*"}' | sort -t '*' -k 9 | cut -d '*' -f1 | cut -d ' ' -f2 | xargs -n 1 dirname`)

VALIDATION_SVN_REV=(`grep '$VDIR' Process_All_Output.sh | grep -v "#" | xargs -n 1 dirname | xargs -n 1 dirname | xargs -n 1 basename | xargs -i svn info {}/FDS_Output_Files | awk '{if($0 != ""){ if(s){s=s"*"$0}else{s=$0}}else{ print s"*";s=""}}END{print s"*"}' | sort -t '*' -k 9 | cut -d '*' -f9 | cut -d ' ' -f4`)

VALIDATION_SVN_DATE=(`grep '$VDIR' Process_All_Output.sh | grep -v "#" | xargs -n 1 dirname | xargs -n 1 dirname | xargs -n 1 basename | xargs -i svn info {}/FDS_Output_Files | awk '{if($0 != ""){ if(s){s=s"*"$0}else{s=$0}}else{ print s"*";s=""}}END{print s"*"}' | sort -t '*' -k 9 | cut -d '*' -f10 | cut -d ' ' -f4`)

# Calculate number of validation sets
NUM_VALIDATION_SETS=${#VALIDATION_SETS[@]}

# Table header
echo "\begin{longtable}[c]{|l|c|c|}" > $OUTPUT_TEX_FILE
echo "\caption[Validation SVN Statistics]{Validation SVN statistics for all data sets}" >> $OUTPUT_TEX_FILE
echo "\label{validation_svn_stats}" >> $OUTPUT_TEX_FILE
echo "\\\\ \hline" >> $OUTPUT_TEX_FILE
echo "Dataset  &  Last Changed SVN Date  &  Last Changed SVN Revision \\\\ \hline \hline" >> $OUTPUT_TEX_FILE
echo "\endfirsthead" >> $OUTPUT_TEX_FILE
echo "\hline" >> $OUTPUT_TEX_FILE
echo "Dataset  &  Last Changed SVN Date  &  Last Changed SVN Revision \\\\ \hline \hline" >> $OUTPUT_TEX_FILE
echo "\endhead" >> $OUTPUT_TEX_FILE

# Loop through data sets and output within table
for i in `seq 0 $(( ${NUM_VALIDATION_SETS} - 1 ))`
do
    echo "${VALIDATION_SETS[i]//_/\_}  &  ${VALIDATION_SVN_DATE[i]}  &  ${VALIDATION_SVN_REV[i]} \\\\ \hline" >> $OUTPUT_TEX_FILE
done

# Table footer
echo "\end{longtable}" >> $OUTPUT_TEX_FILE

cd $CURRENT_DIR

