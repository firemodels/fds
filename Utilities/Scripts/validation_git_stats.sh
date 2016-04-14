#!/bin/bash

# This script outputs a LaTeX file with a table of the FDS validation
# sets and their corresponding GIT information (i.e., when the FDS
# output files were last commited to the repository). This table
# is then included in the FDS Validation Guide.
MAKEGITENTRY(){
DIR=$1
gitrevisions=/tmp/gitrevisions.$$
cat $FDSSMV/Validation/$DIR/FDS_Output_Files/*git.txt 2> /dev/null | sort -u > $gitrevisions
gitrev=`head -1 $gitrevisions`
if [ "$gitrev" != "" ] ; then
  gitrevshort=`echo $gitrev | awk -F - '{print $4}' | sed 's/^.\{1\}//'`
  gitdate=`git show -s --format=%aD $gitrevshort | head -1 | awk '{print $3,$2",",$4}'`
  gitdate2=`git show -s --format=%at $gitrevshort | head -1 | awk '{print $1}'`
  echo "${DIR//_/\_}  & $gitdate & $gitrev & $gitdate2 \\\\ \hline"
fi
rm $gitrevisions
}


CURRENT_DIR=`pwd`
if [ "$FDSSMV" == "" ] ; then
   FDSSMV=~/FDS-SMVgitclean
fi
export FDSSMV

while getopts 'r:' OPTION
do
case $OPTION  in
  r)
   FDSSMV="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

cd $FDSSMV/Utilities/Scripts

# Name and location of output .tex file with validation GIT statistics
OUTPUT_TEX_FILE=$FDSSMV/Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/ScatterPlots/validation_git_stats.tex

# Table header
echo "\begin{longtable}[c]{|l|c|c|}" > $OUTPUT_TEX_FILE
echo "\caption[Validation Git Statistics]{Validation Git statistics for all data sets}" >> $OUTPUT_TEX_FILE
echo "\label{validation_git_stats}" >> $OUTPUT_TEX_FILE
echo "\\\\ \hline" >> $OUTPUT_TEX_FILE
echo "Dataset  &  FDS Revision Date  &  FDS Revision String\\\\ \hline \hline" >> $OUTPUT_TEX_FILE
echo "\endfirsthead" >> $OUTPUT_TEX_FILE
echo "\hline" >> $OUTPUT_TEX_FILE
echo "Dataset  &  FDS Revision Date  &  FDS Revision String\\\\ \hline \hline" >> $OUTPUT_TEX_FILE
echo "\endhead" >> $OUTPUT_TEX_FILE

# Table body
maketable=$FDSSMV/Validation/Process_All_Output.sh
CASELIST=/tmp/temp.out.$$
TABLE_ENTRIES=/tmp/temp2.out.$$
grep VDIR $maketable | awk 'BEGIN { FS = "/" } ; { print $2 }' > $CASELIST
while read p; do
  MAKEGITENTRY   $p  >> $TABLE_ENTRIES
done <$CASELIST
	cat $TABLE_ENTRIES | sort -n -t '&' -k 4 | awk -F "&" '{ print $1 "&" $2 "&" $3 "\\\\ \\hline"}' >> $OUTPUT_TEX_FILE
rm $CASELIST $TABLE_ENTRIES

# Table footer
echo "\end{longtable}" >> $OUTPUT_TEX_FILE

cd $CURRENT_DIR

