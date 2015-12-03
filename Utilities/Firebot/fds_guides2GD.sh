#!/bin/bash
FROMDIR=$1
MANDIR=$2

GDRIVE=~/bin/gdrive
CURDIR=`pwd`
# directory containing guides on google drive : FDS-SMV Newest Manuals
PARENT_ID=0B_wB1pJL2bFQUlJwMmNfaHlqME0

UPLOADGUIDE ()
{
  cd $FROMDIR
  FILE=$1
  FILEnew=${FILE}_new.pdf
  cp $FILE.pdf $FILEnew
  $GDRIVE list  | grep $FILEnew | awk '{ system("~/bin/gdrive delete -i " $1)} '
  $GDRIVE upload -p $PARENT_ID -f $FILEnew
}
UPLOADFIGURES ()
{
  cd $MANDIR/SCRIPT_FIGURES
  DIRECTORY=$1
  tarfile=${DIRECTORY}_figures.tar
  tar cvf ../$tarfile .
  cd ..
  gzip $tarfile
  $GDRIVE list  | grep $tarfile.gz | awk '{ system("~/bin/gdrive delete -i " $1)} '
  $GDRIVE upload -p $PARENT_ID -f $tarfile.gz
}

if [ -e $GDRIVE ] ; then
  UPLOADGUIDE FDS_Config_Management_Plan
  UPLOADGUIDE FDS_Technical_Reference_Guide
  UPLOADGUIDE FDS_User_Guide
  UPLOADGUIDE FDS_Validation_Guide
  UPLOADGUIDE FDS_Verification_Guide
  UPLOADFIGURES FDS_Technical_Reference_Guide
  UPLOADFIGURES FDS_User_Guide
  UPLOADFIGURES FDS_Validation_Guide
  UPLOADFIGURES FDS_Verification_Guide
  cd $CURDIR
fi
