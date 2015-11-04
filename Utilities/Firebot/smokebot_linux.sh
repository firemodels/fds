#!/bin/bash

# Smokebot
# This script is derived from Kris Overholt's firebot script. 
# It tests smokeview by running the smokeview verification suite

#  ===================
#  = Input variables =
#  ===================

# define run directories
SMOKEBOT_RUNDIR=`pwd`
OUTPUT_DIR="$SMOKEBOT_RUNDIR/output"
HISTORY_DIR="$SMOKEBOT_RUNDIR/history"
TIME_LOG=$OUTPUT_DIR/timings
ERROR_LOG=$OUTPUT_DIR/errors
WARNING_LOG=$OUTPUT_DIR/warnings
GUIDE_DIR=$OUTPUT_DIR/guides
STAGE_STATUS=$OUTPUT_DIR/stage_status
NEWGUIDE_DIR=$OUTPUT_DIR/Newest_Guides

# define repo names (default)
fdsroot=~/FDS-SMVgitclean
cfastroot=~/cfastgitclean

SMOKEBOT_QUEUE=smokebot
MAKEMOVIES=
RUNAUTO=
BRANCH=
RUNDEBUG="1"
OPENMP=
RUN_OPENMP=
TESTFLAG=
CLEANREPO=0
UPDATEREPO=0
SSH=
MAILTO=
UPLOADRESULTS=

WEBHOSTNAME=blaze.nist.gov
if [ "$SMOKEBOT_HOSTNAME" != "" ] ; then
WEBHOSTNAME=$SMOKEBOT_HOSTNAME
fi

if [[ "$IFORT_COMPILER" != "" ]] ; then
  source $IFORT_COMPILER/bin/compilervars.sh intel64
fi 
notfound=`icc -help 2>&1 | tail -1 | grep "not found" | wc -l`
if [ "$notfound" == "1" ] ; then
  export haveCC="0"
  USEINSTALL="-i"
  USEINSTALL2="-u"
else
  export haveCC="1"
  USEINSTALL=
  USEINSTALL2=
fi

while getopts 'ab:C:cm:Mo:q:r:sS:tuU' OPTION
do
case $OPTION in
  a)
   RUNAUTO="y"
   ;;
  b)
   BRANCH="$OPTARG"
   ;;
  C)
   cfastroot="$OPTARG"
   ;;
  c)
   CLEANREPO=1
   ;;
  m)
   MAILTO="$OPTARG"
   ;;
  M)
   MAKEMOVIES="1"
   ;;
  o)
   nthreads="$OPTARG"
   OPENMP=openmp_
   RUN_OPENMP="-o $nthreads"
   ;;
  q)
   SMOKEBOT_QUEUE="$OPTARG"
   ;;
  r)
   fdsroot="$OPTARG"
   ;;
  s)
   RUNDEBUG="0"
   ;;
  S)
   SSH="ssh $OPTARG "
   ;;
  t)
   TESTFLAG="-t"
   ;;
  U)
   UPLOADRESULTS=1
   ;;
  u)
   UPDATEREPO=1
   ;;
esac
done
shift $(($OPTIND-1))

DB=_db
IB=
if [ "$FDSNETWORK" == "infiniband" ] ; then
IB=ib
fi

platform="linux"
platform2="Linux"
if [ "`uname`" == "Darwin" ]
then
  platform="osx"
  platform2="OSX"
fi
export platform

cd

export fdsroot
export cfastroot

export SMV_Summary="$fdsroot/Manuals/SMV_Summary"
WEBFROMDIR="$fdsroot/Manuals/SMV_Summary"
WEBTODIR=/var/www/html/VV/SMV2

SMV_VG_GUIDE=$fdsroot/Manuals/SMV_Verification_Guide/SMV_Verification_Guide.pdf
SMV_UG_GUIDE=$fdsroot/Manuals/SMV_User_Guide/SMV_User_Guide.pdf
GEOM_NOTES=$fdsroot/Manuals/FDS_User_Guide/geom_notes.pdf
UploadGuides=$fdsroot/Utilities/Firebot/smv_guides2GD.sh

THIS_FDS_AUTHOR=
THIS_FDS_FAILED=0
FDS_STATUS_FILE=$fdsroot/FDS_status
LAST_FDS_FAILED=0
if [ -e $FDS_STATUS_FILE ] ; then
  LAST_FDS_FAILED=`cat $FDS_STATUS_FILE`
fi

# Load mailing list for status report
source $SMOKEBOT_RUNDIR/firebot_email_list.sh

mailTo=$mailToSMV
if [[ "$LAST_FDS_FAILED" == "1" ]] ; then
  mailTo=$mailToFDS
fi
if [[ "$MAILTO" != "" ]]; then
  mailTo=$MAILTO
fi

JOBPREFIX=SB_

#  =============================================
#  = Smokebot timing and notification mechanism =
#  =============================================

# This routine checks the elapsed time of Smokebot.
# If Smokebot runs more than 12 hours, an email notification is sent.
# This is a notification only and does not terminate Smokebot.
# This check runs during Stages 3 and 5.

# Start firebot timer
START_TIME=$(date +%s)

# Set time limit (43,200 seconds = 12 hours)
TIME_LIMIT=43200
TIME_LIMIT_EMAIL_NOTIFICATION="unsent"

run_auto()
{
  GIT_STATUSDIR=~/.smokebot
  SMV_SOURCE=$fdsroot/SMV/source
  GIT_SMVFILE=$GIT_STATUSDIR/smv_revision
  GIT_SMVLOG=$GIT_STATUSDIR/smv_log

  FDS_SOURCE=$fdsroot/FDS_Source
  GIT_FDSFILE=$GIT_STATUSDIR/fds_revision
  GIT_FDSLOG=$GIT_STATUSDIR/FDS_log

  MESSAGE_FILE=$GIT_STATUSDIR/message

  MKDIR $GIT_STATUSDIR
# remove untracked files, revert repo files, update to latest revision
  cd $fdsroot

  CURRENT_BRANCH=`git rev-parse --abbrev-ref HEAD`
  if [[ "$BRANCH" != "" ]] ; then
    if [[ `git branch | grep $BRANCH` == "" ]] ; then 
       echo "Error: the branch $BRANCH does not exist. Terminating script."
       exit
    fi
    if [[ "$BRANCH" != "$CURRENT_BRANCH" ]] ; then
       echo Checking out branch $BRANCH.
       git checkout $BRANCH
    fi
  else
     BRANCH=$CURRENT_BRANCH
  fi
  if [[ "$UPDATE" == "1" ]] ; then
    echo Update the branch $BRANCH.
    git pull 
  fi

# get info for smokeview
  cd $SMV_SOURCE
  THIS_SMVREVISION=`git log --abbrev-commit . | head -1 | awk '{print $2}'`
  THIS_SMVAUTHOR=`git log . | head -2 | tail -1 | awk '{print $2}'`
  LAST_SMVREVISION=`cat $GIT_SMVFILE`
  git log . | head -5 | tail -1 > $GIT_SMVLOG

# get info for FDS
  cd $FDS_SOURCE
  THIS_FDSREVISION=`git log --abbrev-commit . | head -1 | awk '{printf $2}'`
  THIS_FDSAUTHOR=`git log . | head -2 | tail -1 | awk '{print $2}'`
  LAST_FDSREVISION=`cat $GIT_FDSFILE`
  git log . | head -5 | tail -1 > $GIT_FDSLOG

  if [[ $THIS_SMVREVISION == $LAST_SMVREVISION && $THIS_FDSREVISION == $LAST_FDSREVISION ]] ; then
    exit
  fi

  rm -f $MESSAGE_FILE
  if [[ $THIS_SMVREVISION != $LAST_SMVREVISION ]] ; then
    echo $THIS_SMVREVISION>$GIT_SMVFILE
    echo -e "smokeview source has changed. $LAST_SMVREVISION->$THIS_SMVREVISION($THIS_SMVAUTHOR)" >> $MESSAGE_FILE
    cat $GIT_SMVLOG >> $MESSAGE_FILE
  fi
  if [[ $THIS_FDSREVISION != $LAST_FDSREVISION ]] ; then
    echo $THIS_FDSREVISION>$GIT_FDSFILE
    echo -e "FDS source has changed. $LAST_FDSREVISION->$THIS_FDSREVISION($THIS_FDSAUTHOR)" >> $MESSAGE_FILE
    cat $GIT_FDSLOG >> $MESSAGE_FILE
  fi
  echo -e "Smokebot run initiated." >> $MESSAGE_FILE
  cat $MESSAGE_FILE | mail -s "smokebot run initiated" $mailTo > /dev/null
}

GET_TIME(){
  echo $(date +"%s")
}

GET_DURATION(){
  time_before=$1
  time_after=$2
  DIFF_TIME=`echo $(($time_after-$time_before))`
  TIME_H=`echo $(($DIFF_TIME / 3600 ))`
  TIME_M=`echo $((($DIFF_TIME % 3600 ) / 60))`
  TIME_S=`echo $(($DIFF_TIME % 60 ))`
  if (( "$DIFF_TIME" >= 3600 )) ; then
    echo "${TIME_H}h ${TIME_M}m ${TIME_S}s"
  else
    if (( "$DIFF_TIME" >= 60 )) ; then
      echo "${TIME_M}m ${TIME_S}s"
    else
      echo "${TIME_S}s"
    fi
  fi
}

MKDIR ()
{
  DIR=$1
  if [ ! -d $DIR ]
  then
    echo Creating directory $DIR
    mkdir $DIR
  fi
}

check_time_limit()
{
   if [ "$TIME_LIMIT_EMAIL_NOTIFICATION" == "sent" ]
   then
      # Continue along
      :
   else
      CURRENT_TIME=$(date +%s)
      ELAPSED_TIME=$(echo "$CURRENT_TIME-$START_TIME"|bc)

      if [ $ELAPSED_TIME -gt $TIME_LIMIT ]
      then
         echo -e "smokebot has been running for more than 12 hours in Stage ${TIME_LIMIT_STAGE}. \n\nPlease ensure that there are no problems. \n\nThis is a notification only and does not terminate smokebot." | mail -s "smokebot Notice: smokebot has been running for more than 12 hours." $mailTo > /dev/null
         TIME_LIMIT_EMAIL_NOTIFICATION="sent"
      fi
   fi
}

#  ========================
#  = Additional functions =
#  ========================

set_files_world_readable()
{
   cd $fdsroot
   chmod -R go+r *
}

clean_smokebot_history()
{
   
   # Clean Smokebot metafiles
   MKDIR $SMOKEBOT_RUNDIR > /dev/null
   cd $SMOKEBOT_RUNDIR
   MKDIR guides > /dev/null
   MKDIR history > /dev/null
   MKDIR output > /dev/null
   rm -rf output/* > /dev/null
   MKDIR $NEWGUIDE_DIR > /dev/null
   chmod 775 $NEWGUIDE_DIR
}

#  ========================
#  ========================
#  = Smokebot Build Stages =
#  ========================
#  ========================

#  ===================================
#  = Stage 0 - External dependencies =
#  ===================================

update_and_compile_cfast()
{
   cd $SMOKEBOT_HOME_DIR

   # Check to see if CFAST repository exists
   if [ -e "$cfastroot" ]
   # If yes, then update the CFAST repository and compile CFAST
   then
      if [ "$CLEANREPO" == "1" ]; then
        echo "Cleaning cfast repo:" > $OUTPUT_DIR/stage0_cfast
        cd $cfastroot
        git clean -dxf > /dev/null
        git add . > /dev/null
        git reset --hard HEAD > /dev/null
      fi

      # Update to latest GIT revision
      if [ "$UPDATEREPO" == "1" ]; then
        echo "Updating cfast repo:" >> $OUTPUT_DIR/stage0_cfast
        git pull >> $OUTPUT_DIR/stage0_cfast 2>&1
      fi
   else
      echo "The cfast repo $cfastroot does not exist"
      echo "Aborting  smokebot"
      exit
   fi
    # Build CFAST
    cd $cfastroot/CFAST/intel_${platform}_64
    rm -f cfast7_${platform}_64
    make --makefile ../makefile clean &> /dev/null
    ./make_cfast.sh >> $OUTPUT_DIR/stage0_cfast 2>&1

   # Check for errors in CFAST compilation
   cd $cfastroot/CFAST/intel_${platform}_64
   if [ -e "cfast7_${platform}_64" ]
   then
      stage0_success=true
   else
      echo "Errors from Stage 0 - CFAST:" >> $ERROR_LOG
      echo "CFAST failed to compile" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage0_cfast >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

}

#  ============================
#  = Stage 1 - GIT operations =
#  ============================

clean_git_repo()
{
   # Check to see if FDS repository exists
   if [ -e "$fdsroot" ]
   then
      if [ "$CLEANREPO" == "1" ]; then
        cd $fdsroot
        git clean -dxf > /dev/null
        git add . > /dev/null
        git reset --hard HEAD > /dev/null
      fi
   else
      echo "The FDS repository $fdsroot does not exist." >> $OUTPUT_DIR/stage1 2>&1
      echo "Aborting smokebot" >> $OUTPUT_DIR/stage1 2>&1
      exit
   fi
}

do_git_checkout()
{
   cd $fdsroot

   CURRENT_BRANCH=`git rev-parse --abbrev-ref HEAD`
   if [[ "$BRANCH" != "" ]] ; then
     if [[ `git branch | grep $BRANCH` == "" ]] ; then 
        echo "Error: the branch $BRANCH does not exist."
        echo "Aborting smokebot"
        exit
     fi
     if [[ "$BRANCH" != "$CURRENT_BRANCH" ]] ; then
        echo "Checking out branch $BRANCH." >> $OUTPUT_DIR/stage1 2>&1
        git checkout $BRANCH
     fi
   else
      BRANCH=$CURRENT_BRANCH
   fi
   if [ "$UPDATEREPO" == "1" ]; then
     echo "Updating branch $BRANCH." >> $OUTPUT_DIR/stage1 2>&1
     git pull >> $OUTPUT_DIR/stage1 2>&1
   fi
   GIT_REVISION=`git describe --long --dirty`
}

check_git_checkout()
{
   cd $fdsroot
   # Check for GIT errors
   stage1_success=true
}

#  ==================================
#  = Stage 2a/b - Compile FDS debug =
#  ==================================

compile_fds_mpi_db()
{
   # Clean and compile mpi FDS debug
   cd $fdsroot/FDS_Compilation/mpi_intel_${platform}_64$IB$DB
   rm -f fds_mpi_intel_${platform}_64$IB$DB
   make --makefile ../makefile clean &> /dev/null
   ./make_fds.sh &> $OUTPUT_DIR/stage2b
}

check_compile_fds_mpi_db()
{
   # Check for errors in FDS debug compilation
   cd $fdsroot/FDS_Compilation/mpi_intel_${platform}_64$IB$DB
   if [ -e "fds_mpi_intel_${platform}_64$IB$DB" ]
   then
      stage2b_success=true
   else
      echo "Errors from Stage 2b - Compile FDS MPI debug:" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage2b >> $ERROR_LOG
      echo "" >> $ERROR_LOG
      THIS_FDS_FAILED=1
   fi

   # Check for compiler warnings/remarks
   if [[ `grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage2b| grep -v 'feupdateenv is not implemented'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 2b warnings:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage2b | grep -v 'feupdateenv is not implemented'>> $WARNING_LOG
      echo "" >> $WARNING_LOG
   # if the executable does not exist then an email has already been sent
      if [ -e "fds_mpi_intel_${platform}_64$IB$DB" ] ; then
        THIS_FDS_FAILED=1
      fi
   fi
}

#  ============================================================
#  = Stage 3a - Run Smokeview verification cases (debug mode) =
#  ============================================================

wait_verification_cases_debug_end()
{
   # Scans qstat and waits for verification cases to end
   if [[ "$SMOKEBOT_QUEUE" == "none" ]]
   then
     while [[ `ps -u $USER -f | fgrep .fds | grep -v grep` != '' ]]; do
        JOBS_REMAINING=`ps -u $USER -f | fgrep .fds | grep -v grep | wc -l`
        echo "Waiting for ${JOBS_REMAINING} verification cases to complete." >> $OUTPUT_DIR/stage3a
        TIME_LIMIT_STAGE="3"
        check_time_limit
        sleep 30
     done
   else
     while [[ `qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX` != '' ]]; do
        JOBS_REMAINING=`qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX | wc -l`
        echo "Waiting for ${JOBS_REMAINING} verification cases to complete." >> $OUTPUT_DIR/stage3a
        TIME_LIMIT_STAGE="3"
        check_time_limit
        sleep 30
     done
   fi
}

run_verification_cases_debug()
{
   #  ======================
   #  = Remove .stop files =
   #  ======================

   # Remove all .stop and .err files from Verification directories (recursively)
   if [ "$CLEANREPO" == "1" ]; then
     cd $fdsroot/Verification
     git clean -dxf > /dev/null
   fi

   #  =====================
   #  = Run all SMV cases =
   #  =====================

   cd $fdsroot/Verification/scripts

   # Submit SMV verification cases and wait for them to start
   echo 'Running SMV verification cases:' >> $OUTPUT_DIR/stage3a 2>&1
   ./Run_SMV_Cases.sh $USEINSTALL2 -m 2 -d -q $SMOKEBOT_QUEUE -j $JOBPREFIX >> $OUTPUT_DIR/stage3a 2>&1

   # Wait for SMV verification cases to end
   wait_verification_cases_debug_end

}

check_verification_cases_debug()
{
   # Scan and report any errors in FDS verification cases
   cd $fdsroot/Verification

   if [[ `grep -rIi 'Run aborted' $OUTPUT_DIR/stage3a` == "" ]] && \
      [[ `grep -rIi 'Segmentation' Visualization/* WUI/* Immersed_Boundary_Method/*` == "" ]] && \
      [[ `grep -rI 'ERROR:' Visualization/* WUI/* Immersed_Boundary_Method/*` == "" ]] && \
      [[ `grep -rIi 'STOP: Numerical' Visualization/* WUI/* Immersed_Boundary_Method/*` == "" ]] && \
      [[ `grep -rIi -A 20 'forrtl' Visualization/* WUI/* Immersed_Boundary_Method/*` == "" ]]
   then
      stage3a_success=true
   else
      grep -rIi 'Run aborted' $OUTPUT_DIR/stage3a > $OUTPUT_DIR/stage3a_errors
      grep -rIi 'Segmentation' Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage3a_errors
      grep -rI 'ERROR:' Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage3a_errors
      grep -rIi 'STOP: Numerical' -rIi Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage3a_errors
      grep -rIi -A 20 'forrtl' Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage3a_errors
      
      echo "Errors from Stage 3a - Run verification cases (debug mode):" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage3a_errors >> $ERROR_LOG
      echo "" >> $ERROR_LOG
      THIS_FDS_FAILED=1
   fi
   if [[ `grep 'Warning' -rI $OUTPUT_DIR/stage3a` == "" ]] 
   then
      no_warnings=true
   else
      echo "Stage 3a warnings:" >> $WARNING_LOG
      grep 'Warning' -rI $OUTPUT_DIR/stage3a >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

#  ====================================
#  = Stage 4a/b - Compile FDS release =
#  ====================================

compile_fds_mpi()
{
   # Clean and compile FDS
   cd $fdsroot/FDS_Compilation/mpi_intel_${platform}_64$IB
   rm -f fds_mpi_intel_${platform}_64$IB
   make --makefile ../makefile clean &> /dev/null
   ./make_fds.sh &> $OUTPUT_DIR/stage4b
}

check_compile_fds_mpi()
{
   # Check for errors in FDS compilation
   cd $fdsroot/FDS_Compilation/mpi_intel_${platform}_64$IB
   if [ -e "fds_mpi_intel_${platform}_64$IB" ]
   then
      stage4b_success=true
   else
      echo "Errors from Stage 4b - Compile FDS release:" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage4b >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   # 'performing multi-file optimizations' and 'generating object file' are part of a normal compile
   if [[ `grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage4b | grep -v 'performing multi-file optimizations' | grep -v 'generating object file'| grep -v 'feupdateenv is not implemented'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 4b warnings:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage4b | grep -v 'performing multi-file optimizations' | grep -v 'generating object file'| grep -v 'feupdateenv is not implemented' >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

#  ======================================
#  = Stage 5pre - Compile SMV utilities =
#  ======================================

compile_smv_utilities()
{
   echo "" > $OUTPUT_DIR/stage5pre
   if [ "$haveCC" == "1" ] ; then
   if [ "$SSH" == "" ] ; then 
   # smokeview libraries
   cd $fdsroot/SMV/Build/LIBS/lib_${platform}_intel_64
   echo 'Building Smokeview libraries:' >> $OUTPUT_DIR/stage5pre 2>&1
   ./makelibs.sh >> $OUTPUT_DIR/stage5pre 2>&1

   # smokezip:
   cd $fdsroot/Utilities/smokezip/intel_${platform}_64
   rm -f *.o smokezip_${platform}_64
   echo 'Compiling smokezip:' >> $OUTPUT_DIR/stage5pre 2>&1
   ./make_zip.sh >> $OUTPUT_DIR/stage5pre 2>&1
   echo "" >> $OUTPUT_DIR/stage5pre 2>&1
   
   # smokediff:
   cd $fdsroot/Utilities/smokediff/intel_${platform}_64
   rm -f *.o smokediff_${platform}_64
   echo 'Compiling smokediff:' >> $OUTPUT_DIR/stage5pre 2>&1
   ./make_diff.sh >> $OUTPUT_DIR/stage5pre 2>&1
   echo "" >> $OUTPUT_DIR/stage5pre 2>&1
   
   # background:
   cd $fdsroot/Utilities/background/intel_${platform}_64
   rm -f *.o background
   echo 'Compiling background:' >> $OUTPUT_DIR/stage5pre 2>&1
   ./make_background.sh >> $OUTPUT_DIR/stage5pre 2>&1
   
  # wind2fds:
   cd $fdsroot/Utilities/wind2fds/intel_${platform}_64
   rm -f *.o wind2fds_${platform}_64
   echo 'Compiling wind2fds:' >> $OUTPUT_DIR/stage5pre 2>&1
   ./make_wind.sh >> $OUTPUT_DIR/stage5pre 2>&1
   echo "" >> $OUTPUT_DIR/stage5pre 2>&1
   else
   $SSH \( \
   cd $fdsroot/SMV/Build/LIBS/lib_${platform}_intel_64 \; \
   echo 'Building Smokeview libraries:' >> $OUTPUT_DIR/stage5pre 2>&1 \; \
   ./makelibs.sh >> $OUTPUT_DIR/stage5pre 2>&1 \; \
   cd $fdsroot/Utilities/smokezip/intel_${platform}_64 \; \
   rm -f *.o smokezip_${platform}_64 \; \
   echo 'Compiling smokezip:' >> $OUTPUT_DIR/stage5pre 2>&1 \; \
   ./make_zip.sh >> $OUTPUT_DIR/stage5pre 2>&1 \; \
   echo "" >> $OUTPUT_DIR/stage5pre 2>&1 \; \
   cd $fdsroot/Utilities/smokediff/intel_${platform}_64 \; \
   rm -f *.o smokediff_${platform}_64 \; \
   echo 'Compiling smokediff:' >> $OUTPUT_DIR/stage5pre 2>&1 \; \
   ./make_diff.sh >> $OUTPUT_DIR/stage5pre 2>&1 \; \
   echo "" >> $OUTPUT_DIR/stage5pre 2>&1 \; \
   cd $fdsroot/Utilities/background/intel_${platform}_64 \; \
   rm -f *.o background \; \
   echo 'Compiling background:' >> $OUTPUT_DIR/stage5pre 2>&1 \; \
   ./make_background.sh >> $OUTPUT_DIR/stage5pre 2>&1 \; \
   cd $fdsroot/Utilities/wind2fds/intel_${platform}_64 \; \
   rm -f *.o wind2fds_${platform}_64 \; \
   echo 'Compiling wind2fds:' >> $OUTPUT_DIR/stage5pre 2>&1 \; \
   ./make_wind.sh >> $OUTPUT_DIR/stage5pre 2>&1 \; \
   echo "" >> $OUTPUT_DIR/stage5pre 2>&1  \)
   fi
   else
   echo "Warning: smokeview and utilities not built - C compiler not available" >> $OUTPUT_DIR/stage5pre 2>&1
   fi
}

is_file_installed()
{
  program=$1
  notfound=`$program -help | tail -1 | grep "not found" | wc -l`
  if [ "$notfound" == "1" ] ; then
    stage5pre_success="0"
    echo "***error: $program not installed" >> $OUTPUT_DIR/stage5pre
  fi
}

check_smv_utilities()
{
   if [ "$haveCC" == "1" ] ; then
     # Check for errors in SMV utilities compilation
     cd $fdsroot
     if [ -e "$fdsroot/Utilities/smokezip/intel_${platform}_64/smokezip_${platform}_64" ]  && \
        [ -e "$fdsroot/Utilities/smokediff/intel_${platform}_64/smokediff_${platform}_64" ]  && \
        [ -e "$fdsroot/Utilities/wind2fds/intel_${platform}_64/wind2fds_${platform}_64" ]  && \
        [ -e "$fdsroot/Utilities/background/intel_${platform}_64/background" ]
     then
        stage5pre_success="1"
     else
        stage5pre_success="0"
        echo "Errors from Stage 5pre - Compile SMV utilities:" >> $ERROR_LOG
        cat $OUTPUT_DIR/stage5pre >> $ERROR_LOG
        echo "" >> $ERROR_LOG
     fi
   else
     stage5pre_success="1"
     is_file_installed smokeview
     is_file_installed smokezip
     is_file_installed smokediff
     is_file_installed wind2fds
     is_file_installed background
     if [ "$stage5pre_success" == "0" ] ; then
        echo "Errors from Stage 5pre - Smokeview and utilities:" >> $ERROR_LOG
        stage5pre_success="1"
        cat $OUTPUT_DIR/stage5pre >> $ERROR_LOG
        echo "" >> $ERROR_LOG
     fi
   fi
}

#  ===================================================
#  = Stage 5 - Run verification cases (release mode) =
#  ===================================================

wait_verification_cases_release_end()
{
   # Scans qstat and waits for verification cases to end
   if [[ "$SMOKEBOT_QUEUE" == "none" ]]
   then
     while [[ `ps -u $USER -f | fgrep .fds | grep -v grep` != '' ]]; do
        JOBS_REMAINING=`ps -u $USER -f | fgrep .fds | grep -v grep | wc -l`
        echo "Waiting for ${JOBS_REMAINING} verification cases to complete." >> $OUTPUT_DIR/stage5
        TIME_LIMIT_STAGE="5"
        check_time_limit
        sleep 60
     done
   else
     while [[ `qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX` != '' ]]; do
        JOBS_REMAINING=`qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX | wc -l`
        echo "Waiting for ${JOBS_REMAINING} verification cases to complete." >> $OUTPUT_DIR/stage5
        TIME_LIMIT_STAGE="5"
        check_time_limit
        sleep 60
     done
   fi
}

run_verification_cases_release()
{
   #  ======================
   #  = Remove .stop files =
   #  ======================

   # Remove all .stop and .err files from Verification directories (recursively)
   if [ "$CLEANREPO" == "1" ]; then
     cd $fdsroot/Verification
     git clean -dxf > /dev/null
   fi

   # Start running all SMV verification cases
   cd $fdsroot/Verification/scripts
   echo 'Running SMV verification cases:' >> $OUTPUT_DIR/stage5 2>&1
   ./Run_SMV_Cases.sh $USEINSTALL2 $RUN_OPENMP -q $SMOKEBOT_QUEUE -j $JOBPREFIX >> $OUTPUT_DIR/stage5 2>&1

   # Wait for all verification cases to end
   wait_verification_cases_release_end
}

check_verification_cases_release()
{
   # Scan and report any errors in FDS verification cases
   cd $fdsroot/Verification

   if [[ `grep -rIi 'Run aborted' $OUTPUT_DIR/stage5` == "" ]] && \
      [[ `grep -rIi 'Segmentation' Visualization/* WUI/* Immersed_Boundary_Method/* ` == "" ]] && \
      [[ `grep -rI 'ERROR:' Visualization/* WUI/* Immersed_Boundary_Method/* ` == "" ]] && \
      [[ `grep -rIi 'STOP: Numerical' Visualization/* WUI/* Immersed_Boundary_Method/* ` == "" ]] && \
      [[ `grep -rIi  -A 20 'forrtl' Visualization/* WUI/* Immersed_Boundary_Method/* ` == "" ]]
   then
      stage5_success=true
   else
      grep -rIi 'Run aborted' $OUTPUT_DIR/stage5 > $OUTPUT_DIR/stage5_errors
      grep -rIi 'Segmentation' Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage5_errors
      grep -rI 'ERROR:' Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage5_errors
      grep -rIi 'STOP: Numerical' Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage5_errors
      grep -rIi -A 20 'forrtl' Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage5_errors

      echo "Errors from Stage 5 - Run verification cases (release mode):" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage5_errors >> $ERROR_LOG
      echo "" >> $ERROR_LOG
      THIS_FDS_FAILED=1
   fi

      
   if [[ `grep 'Warning' -rI $OUTPUT_DIR/stage5` == "" ]] 
   then
      no_warnings=true
   else
      echo "Stage 5 warnings:" >> $WARNING_LOG
      grep 'Warning' -rI $OUTPUT_DIR/stage5 >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

#  ================================
#  = Stage 6a - Compile SMV debug =
#  ================================

compile_smv_db()
{
   if [ "$haveCC" == "1" ] ; then
   if [ "$SSH" == "" ] ; then
   # Clean and compile SMV debug
   cd $fdsroot/SMV/Build/intel_${platform}_64
   rm -f smokeview_${platform}_64_db
   ./make_smv_db.sh &> $OUTPUT_DIR/stage6a
   else
   $SSH \(
   cd $fdsroot/SMV/Build/intel_${platform}_64 \; \
   rm -f smokeview_${platform}_64_db \; \
   ./make_smv_db.sh &> $OUTPUT_DIR/stage6a \)
   fi
   fi
}

check_compile_smv_db()
{
   if [ "$haveCC" == "1" ] ; then
   # Check for errors in SMV debug compilation
   cd $fdsroot/SMV/Build/intel_${platform}_64
   if [ -e "smokeview_${platform}_64_db" ]
   then
      stage6a_success=true
   else
      echo "Errors from Stage 6a - Compile SMV debug:" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage6a >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   # grep -v 'feupdateenv ...' ignores a known FDS MPI compiler warning (http://software.intel.com/en-us/forums/showthread.php?t=62806)
   if [[ `grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage6a | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 6a warnings:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage6a | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked' >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
   fi
}

#  =============================================
#  = Stage 6b - Make SMV pictures (debug mode) =
#  =============================================

make_smv_pictures_db()
{
   # Run Make SMV Pictures script (debug mode)
   if [ "$SSH" == "" ]; then
   cd $fdsroot/Verification/scripts
   ./Make_SMV_Pictures.sh $USEINSTALL -d 2>&1 | grep -v FreeFontPath &> $OUTPUT_DIR/stage6b
   else
   $SSH \( cd $fdsroot/Verification/scripts \; \
   ./Make_SMV_Pictures.sh $USEINSTALL -d 2>&1 \| grep -v FreeFontPath &> $OUTPUT_DIR/stage6b \)
   fi
}

check_smv_pictures_db()
{
   # Scan and report any errors in make SMV pictures process
   cd $SMOKEBOT_RUNDIR
   if [[ `grep -I -E "Segmentation|Error" $OUTPUT_DIR/stage6b` == "" ]]
   then
      stage6b_success=true
   else
      cp $OUTPUT_DIR/stage6b $OUTPUT_DIR/stage6b_errors

      echo "Errors from Stage 6b - Make SMV pictures (debug mode):" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage6b_errors >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Scan for and report any warnings in make SMV pictures process
   cd $SMOKEBOT_RUNDIR
   if [[ `grep -I -E "Warning" $OUTPUT_DIR/stage6b` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 6b - Make SMV pictures (debug mode):" >> $WARNING_LOG
      grep -I -E "Warning" $OUTPUT_DIR/stage6b >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi

}

#  ==================================
#  = Stage 6c - Compile SMV release =
#  ==================================

compile_smv()
{
   if [ "$haveCC" == "1" ] ; then
   if [ "$SSH" == "" ] ; then
   # Clean and compile SMV
   cd $fdsroot/SMV/Build/intel_${platform}_64
   rm -f smokeview_${platform}_64
   ./make_smv.sh $TESTFLAG &> $OUTPUT_DIR/stage6c
   else
   $SSH \( \
   cd $fdsroot/SMV/Build/intel_${platform}_64 \; \
   rm -f smokeview_${platform}_64 \; \
   ./make_smv.sh $TESTFLAG &> $OUTPUT_DIR/stage6c \)
   fi
   fi
}

check_compile_smv()
{
   if [ "$haveCC" == "1" ] ; then
   # Check for errors in SMV release compilation
   cd $fdsroot/SMV/Build/intel_${platform}_64
   if [ -e "smokeview_${platform}_64" ]
   then
      stage6c_success=true
   else
      echo "Errors from Stage 6c - Compile SMV release:" >> $ERROR_LOG
      echo "The program smokeview_${platform}_64 does not exist."
      cat $OUTPUT_DIR/stage6c >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   # grep -v 'feupdateenv ...' ignores a known FDS MPI compiler warning (http://software.intel.com/en-us/forums/showthread.php?t=62806)
   if [[ `grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage6c | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 6c warnings:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage6c | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked' >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
   fi
}

#  ===============================================
#  = Stage 6d - Make SMV pictures (release mode) =
#  ===============================================

make_smv_pictures()
{
   # Run Make SMV Pictures script (release mode)
   if [ "$SSH" == "" ]; then
   cd $fdsroot/Verification/scripts
   ./Make_SMV_Pictures.sh $TESTFLAG $USEINSTALL 2>&1 | grep -v FreeFontPath &> $OUTPUT_DIR/stage6d
   else
   $SSH \( cd $fdsroot/Verification/scripts \; \
   ./Make_SMV_Pictures.sh $TESTFLAG $USEINSTALL \| grep -v FreeFontPath &> $OUTPUT_DIR/stage6d \)
   fi
}

check_smv_pictures()
{
   # Scan and report any errors in make SMV pictures process
   cd $SMOKEBOT_RUNDIR
   if [[ `grep -I -E "Segmentation|Error" $OUTPUT_DIR/stage6d` == "" ]]
   then
      stage6d_success=true
   else
      cp $OUTPUT_DIR/stage6d  $OUTPUT_DIR/stage6d_errors

      echo "Errors from Stage 6d - Make SMV pictures (release mode):" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage6d >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi
}

#  ===============================================
#  = Stage 6e - Make SMV movies (release mode) =
#  ===============================================

make_smv_movies()
{
   cd $fdsroot/Verification
   scripts/Make_SMV_Movies.sh 2>&1  &> $OUTPUT_DIR/stage6e
}

check_smv_movies()
{
   cd $SMOKEBOT_RUNDIR
   if [[ `grep -I -E "Segmentation|Error" $OUTPUT_DIR/stage6e` == "" ]]
   then
      stage6e_success=true
   else
      cp $OUTPUT_DIR/stage6e  $OUTPUT_DIR/stage6e_errors

      echo "Errors from Stage 6e - Make SMV movies " >> $ERROR_LOG
      cat $OUTPUT_DIR/stage6e >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Scan for and report any warnings in make SMV pictures process
   cd $SMOKEBOT_RUNDIR
   if [[ `grep -I -E "Warning" $OUTPUT_DIR/stage6e` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 6e - Make SMV movies (release mode):" >> $WARNING_LOG
      grep -I -E "Warning" $OUTPUT_DIR/stage6e >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
   if [ "$UPLOADRESULTS" == "1" ]; then
     if [ -d "$WEBTODIR" ]; then
       if [ -d "$WEBFROMDIR" ]; then 
         CURDIR=`pwd`
         cd $WEBTODIR
         rm -rf *
         cd $WEBFROM
         cp -r * $WEBTODIR/.
         cd $CURDIR
       fi
     fi
   fi

}

#  ======================================
#  = Stage 7 - FDS run time statistics =
#  ======================================

generate_timing_stats()
{
   cd $fdsroot/Verification/scripts/
   export QFDS="$fdsroot/Verification/scripts/copyout.sh"
   export RUNCFAST="$fdsroot/Verification/scripts/copyout.sh"
   export RUNTFDS="$fdsroot/Verification/scripts/copyout.sh"

   cd $fdsroot/Verification
   scripts/SMV_Cases.sh
   scripts/SMV_geom_Cases.sh

   cd $fdsroot/Utilities/Scripts
   ./fds_timing_stats.sh smokebot
}

archive_timing_stats()
{
   cd $fdsroot/Utilities/Scripts
   cp fds_timing_stats.csv "$HISTORY_DIR/${GIT_REVISION}_timing.csv"
}

#  ===================================
#  = Stage 8 - Build smokview guides =
#  ===================================

check_guide()
{
   stage=$1
   directory=$2
   document=$3
   label=$4

   # Scan and report any errors in build process for guides
   SMOKEBOT_MANDIR=/var/www/html/smokebot/manuals/
   cd $SMOKEBOT_RUNDIR
   if [[ `grep "! LaTeX Error:" -I $stage` == "" ]]
   then
      if [ "$UPLOADRESULTS" == "1" ]; then
      if [ -d $SMOKEBOT_MANDIR ] ; then
        cp $directory/$document $SMOKEBOT_MANDIR/.
      fi
      fi
      if [ -d $SMV_Summary/manuals ] ; then
        cp $directory/$document $SMV_Summary/manuals/.
      fi
      cp $directory/$document $NEWGUIDE_DIR/.
      chmod 664 $NEWGUIDE_DIR/$document
   else
      echo "Errors from Stage 8 - Build FDS-SMV Guides:" >> $ERROR_LOG
      echo $3 >> $ERROR_LOG
      grep "! LaTeX Error:" -I $1 >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for LaTeX warnings (undefined references or duplicate labels)
   if [[ `grep -E "undefined|multiply defined|multiply-defined" -I ${stage}` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 8 warnings:" >> $WARNING_LOG
      echo $label >> $WARNING_LOG
      grep -E "undefined|multiply defined|multiply-defined" -I $stage >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

make_guide()
{
   document=$1
   directory=$2
   label=$3

   cd $directory
  
   ./make_guide.sh &> $OUTPUT_DIR/stage8_$document

   # Check guide for completion and copy to website if successful
   check_guide $OUTPUT_DIR/stage8_$document $directory $document.pdf $label
}

#  =====================================================
#  = Build status reporting - email and save functions =
#  =====================================================

save_build_status()
{
   cd $SMOKEBOT_RUNDIR
   # Save status outcome of build to a text file
   if [[ -e $WARNING_LOG && -e $ERROR_LOG ]]
   then
     cat "" >> $ERROR_LOG
     cat $WARNING_LOG >> $ERROR_LOG
     echo "Build failure and warnings for Version: ${GIT_REVISION}, Branch: $BRANCH." > "$HISTORY_DIR/${GIT_REVISION}.txt"
     cat $ERROR_LOG > "$HISTORY_DIR/${GIT_REVISION}_errors.txt"
     touch output/status_errors_and_warnings

   # Check for errors only
   elif [ -e $ERROR_LOG ]
   then
      echo "Build failure for Version: ${GIT_REVISION}, Branch: $BRANCH." > "$HISTORY_DIR/${GIT_REVISION}.txt"
      cat $ERROR_LOG > "$HISTORY_DIR/${GIT_REVISION}_errors.txt"
      touch output/status_errors

   # Check for warnings only
   elif [ -e $WARNING_LOG ]
   then
      echo "Version: ${GIT_REVISION}, Branch: $BRANCH has warnings." > "$HISTORY_DIR/${GIT_REVISION}.txt"
      cat $WARNING_LOG > "$HISTORY_DIR/${GIT_REVISION}_warnings.txt"
      touch output/status_warnings

   # No errors or warnings
   else
      echo "Build success! Version: ${GIT_REVISION}, Branch: $BRANCH passed all build tests." > "$HISTORY_DIR/${GIT_REVISION}.txt"
      touch output/status_success
   fi
}

email_build_status()
{
   if [[ "$THIS_FDS_FAILED" == "1" ]] ; then
     mailTo="$mailToFDS"
   fi
   echo $THIS_FDS_FAILED>$FDS_STATUS_FILE
   stop_time=`date`
   echo "----------------------------------------------" > $TIME_LOG
   echo ".         host: $hostname " >> $TIME_LOG
   echo ".        start: $start_time " >> $TIME_LOG
   echo ".         stop: $stop_time " >> $TIME_LOG
   echo ".    run cases: $DIFF_RUNCASES" >> $TIME_LOG
   echo ".make pictures: $DIFF_MAKEPICTURES" >> $TIME_LOG
if [ "$MAKEMOVIES" == "1" ]; then
   echo ".  make movies: $DIFF_MAKEMOVIES" >> $TIME_LOG
fi
   echo ".        total: $DIFF_SCRIPT_TIME" >> $TIME_LOG
if [ "$RUNAUTO" != "" ]; then
   echo ".FDS revisions: old: $LAST_FDSREVISION new: $THIS_FDSREVISION" >> $TIME_LOG
   echo ".SMV revisions: old: $LAST_SMVREVISION new: $THIS_SMVREVISION" >> $TIME_LOG
fi
  if [[ $THIS_SMVREVISION != $LAST_SMVREVISION ]] ; then
    cat $GIT_SMVLOG >> $TIME_LOG
  fi
  if [[ $THIS_FDSREVISION != $LAST_FDSREVISION ]] ; then
    cat $GIT_FDSLOG >> $TIME_LOG
  fi
   echo "----------------------------------------------" >> $TIME_LOG
   cd $SMOKEBOT_RUNDIR
   # Check for warnings and errors
   echo "Nightly Manuals (private): http://$WEBHOSTNAME/VV/SMV2" >> $TIME_LOG
   echo "Nightly Manuals  (public):  http://goo.gl/n1Q3WH" >> $TIME_LOG
   echo "-------------------------------" >> $TIME_LOG
   if [[ -e $WARNING_LOG && -e $ERROR_LOG ]]
   then
     # Send email with failure message and warnings, body of email contains appropriate log file
     cat $ERROR_LOG $TIME_LOG | mail -s "smokebot build failure and warnings on ${hostname}. Version: ${GIT_REVISION}, Branch: $BRANCH." $mailTo > /dev/null

   # Check for errors only
   elif [ -e $ERROR_LOG ]
   then
      # Send email with failure message, body of email contains error log file
      cat $ERROR_LOG $TIME_LOG | mail -s "smokebot build failure on ${hostname}. Version: ${GIT_REVISION}, Branch: $BRANCH." $mailTo > /dev/null

   # Check for warnings only
   elif [ -e $WARNING_LOG ]
   then
     # Send email with success message, include warnings
     cat $WARNING_LOG $TIME_LOG | mail -s "smokebot build success with warnings on ${hostname}. Version: ${GIT_REVISION}, Branch: $BRANCH." $mailTo > /dev/null

   # No errors or warnings
   else
# upload guides to a google drive directory
      if [ "$UPLOADRESULTS" == "1" ];then
        cd $SMOKEBOT_RUNDIR
        $UploadGuides $NEWGUIDE_DIR > /dev/null
      fi

      # Send success message with links to nightly manuals
      cat $TIME_LOG | mail -s "smokebot build success on ${hostname}! Version: ${GIT_REVISION}, Branch: $BRANCH." $mailTo > /dev/null
   fi
}

# if -a option is invoked, only proceed running smokebot if the
# smokeview or FDS source has changed

if [[ $RUNAUTO == "y" ]] ; then
  run_auto
fi

#  ============================
#  = Primary script execution =
#  ============================

SCRIPT_TIME_beg=`GET_TIME`
PRELIM_beg=`GET_TIME`
echo "" > $STAGE_STATUS
hostname=`hostname`
start_time=`date`
clean_smokebot_history

### Stage 0 ###
update_and_compile_cfast

### Stage 1 ###
clean_git_repo
do_git_checkout
check_git_checkout
PRELIM_end=`GET_TIME`
DIFF_PRELIM=`GET_DURATION $PRELIM_beg $PRELIM_end`
echo "Preliminary: $DIFF_PRELIM" >> $STAGE_STATUS

### Stage 2b ###
BUILDFDS_beg=`GET_TIME`
compile_fds_mpi_db
check_compile_fds_mpi_db

### Stage 4b ###
stage4_beg=`GET_TIME`
if [[ $stage2b_success ]] ; then
   compile_fds_mpi
   check_compile_fds_mpi
fi
BUILDFDS_end=`GET_TIME`
DIFF_BUILDFDS=`GET_DURATION $BUILDFDS_beg $BUILDFDS_end`
echo "Build FDS: $DIFF_BUILDFDS" >> $STAGE_STATUS

### Stage 5pre ###
SMVUTILSpre_beg=`GET_TIME`
compile_smv_utilities
check_smv_utilities
SMVUTILSpre_end=`GET_TIME`
DIFF_SMVUTILSpre=`GET_DURATION $SMVUTILSpre_beg $SMVUTILSpre_end`
echo "Build SMV Utilities: $DIFF_SMVUTILSpre" >> $STAGE_STATUS

### Stage 3 ###
RUNCASES_beg=`GET_TIME`
if [[ $stage2b_success && "$RUNDEBUG" == "1" ]] ; then
   run_verification_cases_debug
   check_verification_cases_debug
fi

### Stage 5 ###
if [[ $stage4b_success ]] ; then
   run_verification_cases_release
   check_verification_cases_release
fi
RUNCASES_end=`GET_TIME`
DIFF_RUNCASES=`GET_DURATION $RUNCASES_beg $RUNCASES_end`
echo "Run cases: $DIFF_RUNCASES" >> $STAGE_STATUS

### Stage 6a ###
BUILDSMV_beg=`GET_TIME`
compile_smv_db
check_compile_smv_db

### Stage 6c ###
compile_smv
check_compile_smv
BUILDSMV_end=`GET_TIME`
DIFF_BUILDSMV=`GET_DURATION $BUILDSMV_beg $BUILDSMV_end`
echo "Build SMV: $DIFF_BUILDSMV" >> $STAGE_STATUS

### Stage 6d ###
MAKEPICTURES_beg=`GET_TIME`
if [[ $stage4b_success && $stage6c_success ]] ; then
  make_smv_pictures
  check_smv_pictures
fi
MAKEPICTURES_end=`GET_TIME`
DIFF_MAKEPICTURES=`GET_DURATION $MAKEPICTURES_beg $MAKEPICTURES_end`
echo "Make pictures: $DIFF_MAKEPICTURES" >> $STAGE_STATUS

### Stage 6e ###
if [ "$MAKEMOVIES" == "1" ]
then
  MAKEMOVIES_beg=`GET_TIME`
 
  make_smv_movies
  check_smv_movies

  MAKEMOVIES_end=`GET_TIME`
  DIFF_MAKEMOVIES=`GET_DURATION $MAKEMOVIES_beg $MAKEMOVIES_end`
  echo "Make movies: $DIFF_MAKEMOVIES" >> $STAGE_STATUS
fi

### Stage 7 ###
if [[ $stage4b_success ]] ; then
  generate_timing_stats
  archive_timing_stats
fi

### Stage 8 ###
MAKEGUIDES_beg=`GET_TIME`
if [[ $stage4b_success && $stage6d_success ]] ; then
#  make_guide geom_notes $fdsroot/Manuals/FDS_User_Guide 'geometry notes'
  make_guide SMV_User_Guide $fdsroot/Manuals/SMV_User_Guide 'SMV User Guide'
  make_guide SMV_Technical_Reference_Guide $fdsroot/Manuals/SMV_Technical_Reference_Guide 'SMV Technical Reference Guide'
  make_guide SMV_Verification_Guide $fdsroot/Manuals/SMV_Verification_Guide 'SMV Verification Guide'
fi
MAKEGUIDES_end=`GET_TIME`
DIFF_MAKEGUIDES=`GET_DURATION $MAKEGUIDES_beg $MAKEGUIDES_end`
echo "Make guides: $DIFF_MAKEGUIDES" >> $STAGE_STATUS

SCRIPT_TIME_end=`GET_TIME`
DIFF_SCRIPT_TIME=`GET_DURATION $SCRIPT_TIME_beg $SCRIPT_TIME_end`
echo "Total time: $DIFF_SCRIPT_TIME" >> $STAGE_STATUS

### Report results ###
set_files_world_readable
save_build_status
email_build_status
