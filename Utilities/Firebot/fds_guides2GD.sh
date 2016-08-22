#!/bin/bash
FROMDIR=$1
MANDIR=$2

GDRIVE=~/bin/gdrive
CURDIR=`pwd`
# directory containing guides on google drive : FDS-SMV Newest Manuals
MANUAL_PARENT_ID=0B_wB1pJL2bFQUlJwMmNfaHlqME0
FIGURES_PARENT_ID=0B-W-dkXwdHWNOGVsZXNzTjdLek0

UPLOADGUIDE ()
{
  cd $FROMDIR
  FILE=$1
  FILEnew=${FILE}.pdf
  $GDRIVE list  | grep $FILEnew | awk '{ system("~/bin/gdrive delete -i " $1)} '
  $GDRIVE upload -p $MANUAL_PARENT_ID -f $FILEnew
}
UPLOADFIGURES ()
{
  DIRECTORY=$1
  FILE=$2
  cd $MANDIR/$DIRECTORY/SCRIPT_FIGURES
  tarfile=${FILE}_figures.tar
  rm -f ../$tarfile
  rm -f ../$tarfile.gz
  tar cvf ../$tarfile . &> /dev/null
  cd ..
  gzip $tarfile
  $GDRIVE list  | grep $tarfile.gz | awk '{ system("~/bin/gdrive delete -i " $1)} '
  $GDRIVE upload -p $FIGURES_PARENT_ID -f $tarfile.gz
}

if [ -e $GDRIVE ] ; then
  UPLOADGUIDE FDS_Config_Management_Plan
  UPLOADGUIDE FDS_Technical_Reference_Guide
  UPLOADGUIDE FDS_User_Guide
  UPLOADGUIDE FDS_Validation_Guide
  UPLOADGUIDE FDS_Verification_Guide
  UPLOADFIGURES FDS_Technical_Reference_Guide FDS_TG
  UPLOADFIGURES FDS_User_Guide FDS_UG
  UPLOADFIGURES FDS_Validation_Guide FDS_VALG
  UPLOADFIGURES FDS_Verification_Guide FDS_VERG
  cd $CURDIR
fi
