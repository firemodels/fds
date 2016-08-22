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
  cd $FROMDIR
  UPLOADGUIDE SMV_User_Guide
  UPLOADGUIDE SMV_Technical_Reference_Guide
  UPLOADGUIDE SMV_Verification_Guide
  UPLOADFIGURES SMV_User_Guide SMV_UG
  UPLOADFIGURES SMV_Verification_Guide SMV_VG
  cd $CURDIR
fi
