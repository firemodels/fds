#!/bin/bash

# This script outputs a LaTeX file with a table of the FDS validation
# sets and their corresponding GIT information (i.e., when the FDS
# output files were last commited to the repository). This table
# is then included in the FDS Validation Guide.
MAKEGITENTRY(){
DIR=$1
gitrevisions=$TEMPDIR/gitrevisions.$$
cat $FIREMODELS_ROOT/out/$DIR/*git.txt 2> /dev/null | sort -u > $gitrevisions
gitrev=`head -1 $gitrevisions`
if [ "$gitrev" != "" ] ; then
  gitrevshort=`echo $gitrev | awk -F - '{print $(NF-1)}' | sed 's/^g//'`
  gitdate=`git show -s --format=%aD $gitrevshort 2> /dev/null | head -1 | awk '{print $3,$2",",$4}'`
  if [ "$gitdate" == "" ]; then
    gitdate="undefined"
    gitdate2=2000000000
    if [ -e ~/FDS-SMV ]; then
      CUR_DIR=`pwd`
      cd ~/FDS-SMV
      gitrevshort=`echo $gitrev | awk -F - '{print $4}' | sed 's/^.\{1\}//'`
      gitdateold=`git show -s --format=%aD $gitrevshort 2> /dev/null | head -1 | awk '{print $3,$2",",$4}'`
      if [ "$gitdateold" != "" ]; then
        gitdate=$gitdateold
        gitdate2=`git show -s --format=%at $gitrevshort | head -1 | awk '{print $1}'`
      fi
      cd $CUR_DIR
    fi
  else
    gitdate2=`git show -s --format=%at $gitrevshort | head -1 | awk '{print $1}'`
  fi
  echo "${DIR//_/\\_}  & $gitdate & $gitrev & $gitdate2 \\\\ \hline"
fi
rm $gitrevisions
}

CURRENT_DIR=`pwd`

#determine repo root 
SCRIPTDIR=`dirname "$(readlink -f "$0")"`
cd $SCRIPTDIR/../../..
FIREMODELS_ROOT=`pwd`
cd $CURDIR

TEMPDIR=$HOME/temp
if [ ! -d $TEMPDIR ]; then
   mkdir $TEMPDIR
fi

export FIREMODELS_ROOT

# the following is kept so older versions of firebot won't crash (ie firebot's that call this script using the -r parameter)
while getopts 'r:' OPTION
do
case $OPTION  in
  r)
   IGNORE="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

cd $FIREMODELS_ROOT/fds/Utilities/Scripts

# Name and location of output .tex file with validation GIT statistics
OUTPUT_TEX_FILE=$FIREMODELS_ROOT/fds/Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/ScatterPlots/validation_git_stats.tex

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
maketable=$FIREMODELS_ROOT/fds/Validation/Process_All_Output.sh
CASELIST=$TEMPDIR/temp.out.$$
TABLE_ENTRIES=$TEMPDIR/temp2.out.$$
grep PROCESS $maketable | awk 'BEGIN { FS = " " } ; { print $2 }' | awk '{if(NR>1)print}'> $CASELIST
while read p; do
  MAKEGITENTRY   $p  >> $TABLE_ENTRIES
done <$CASELIST
	cat $TABLE_ENTRIES | sort -n -t '&' -k 4 | awk -F "&" '{ print $1 "&" $2 "&" $3 "\\\\ \\hline"}' >> $OUTPUT_TEX_FILE
rm $CASELIST $TABLE_ENTRIES


# Table footer
echo "\end{longtable}" >> $OUTPUT_TEX_FILE

cd $CURRENT_DIR

