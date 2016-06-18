#!/bin/bash

# Smokebot
# This script is derived from Kris Overholt's firebot script. 
# It tests smokeview by running the smokeview verification suite

#  ===================
#  = Input variables =
#  ===================

size=_64
# define run directories
SMOKEBOT_RUNDIR=`pwd`
OUTPUT_DIR="$SMOKEBOT_RUNDIR/output"
HISTORY_DIR="$HOME/.smokebot/history"
TIME_LOG=$OUTPUT_DIR/timings
ERROR_LOG=$OUTPUT_DIR/errors
WARNING_LOG=$OUTPUT_DIR/warnings
GUIDE_DIR=$OUTPUT_DIR/guides
STAGE_STATUS=$OUTPUT_DIR/stage_status
NEWGUIDE_DIR=$OUTPUT_DIR/Newest_Guides
web_DIR=
WEB_URL=
SMOKEBOT_LITE=

# define repo names (default)
fdsrepo=~/FDS-SMVgitclean
cfastrepo=~/cfastgitclean

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
COMPILER=intel

while getopts 'aAb:C:cI:Lm:Mo:q:r:sS:tuUw:W:' OPTION
do
case $OPTION in
  a)
   RUNAUTO="y"
   ;;
  A)
   RUNAUTO="Y"
   ;;
  b)
   BRANCH="$OPTARG"
   ;;
  C)
   cfastrepo="$OPTARG"
   ;;
  c)
   CLEANREPO=1
   ;;
  I)
   COMPILER="$OPTARG"
   ;;
  L)
   SMOKEBOT_LITE=1
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
   fdsrepo="$OPTARG"
   ;;
  s)
   RUNDEBUG="0"
   ;;
  S)
   SSH="$OPTARG"
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
  w)
   web_DIR="$OPTARG"
   ;;
  W)
   WEB_URL="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

# if one of WEB_URL or web_DIR exist then both should exist
# if web_DIR exists then it must be writable

if [ "$WEB_URL" == "" ]; then
  web_DIR=
fi
if [ "$web_DIR" == "" ]; then
  WEB_URL=
else
  if [ -d $web_DIR ]; then
    testfile=$web_DIR/test.$$
    touch $testfile >& /dev/null
    if [ -e $testfile ]; then
      rm $testfile
    else
      web_DIR=
      WEB_URL=
    fi
  else
    web_DIR=
    WEB_URL=
  fi
fi

if [ "$COMPILER" == "intel" ]; then
if [[ "$IFORT_COMPILER" != "" ]] ; then
  source $IFORT_COMPILER/bin/compilervars.sh intel64
fi 
notfound=`icc -help 2>&1 | tail -1 | grep "not found" | wc -l`
else
notfound=`gcc -help 2>&1 | tail -1 | grep "not found" | wc -l`
fi
if [ "$notfound" == "1" ] ; then
  export haveCC="0"
  USEINSTALL="-i"
  USEINSTALL2="-u"
else
  export haveCC="1"
  USEINSTALL=
  USEINSTALL2=
fi

if [ "$SSH" != "" ]; then
  sshok=$(ssh -o BatchMode=yes -o ConnectTimeout=5 $SSH echo ok 2>/dev/null)
  if [ "$sshok" != "ok" ]; then
    echo unable to make an ssh connection to $SSH
    echo smokebot aborted
    exit
  fi
  SSH="ssh $SSH "
fi

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

echo ""
echo "Preliminaries:"
echo "  running in: $SMOKEBOT_RUNDIR"
echo "FDS-SMV repo: $fdsrepo"
echo "  cfast repo: $cfastrepo"
if [ ! "$web_DIR" == "" ]; then
echo "     web dir: $web_DIR"
fi
if [ ! "$WEB_URL" == "" ]; then
echo "         URL: $WEB_URL"
fi
echo ""

cd

export fdsrepo
export cfastrepo

export SMV_SUMMARY="$fdsrepo/Manuals/SMV_Summary"
WEBFROMDIR="$fdsrepo/Manuals/SMV_Summary"

SMV_VG_GUIDE=$fdsrepo/Manuals/SMV_Verification_Guide/SMV_Verification_Guide.pdf
SMV_UG_GUIDE=$fdsrepo/Manuals/SMV_User_Guide/SMV_User_Guide.pdf
GEOM_NOTES=$fdsrepo/Manuals/FDS_User_Guide/geom_notes.pdf
UploadGuides=$fdsrepo/Utilities/Firebot/smv_guides2GD.sh

THIS_FDS_AUTHOR=
THIS_FDS_FAILED=0
THIS_CFAST_FAILED=0
FDS_STATUS_FILE=$fdsrepo/FDS_status
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
  option=$1
  GIT_STATUSDIR=~/.smokebot
  SMV_SOURCE=$fdsrepo/SMV/source
  TRIGGER_DIR=$fdsrepo/SMV/source/scripts
  GIT_SMV_FILE=$GIT_STATUSDIR/smv_revision
  GIT_SMV_LOG=$GIT_STATUSDIR/smv_log
  
  QUICKTRIGGER=$TRIGGER_DIR/smokeview/smokebot_quicktrigger.txt
  GIT_QT_FILE=$GIT_STATUSDIR/quicktrigger_revision

  TRIGGER=$TRIGGER_DIR/smokeview/smokebot_trigger.txt
  GIT_T_FILE=$GIT_STATUSDIR/trigger_revision

  FDS_SOURCE=$fdsrepo/FDS_Source
  GIT_FDS_FILE=$GIT_STATUSDIR/fds_revision
  GIT_FDS_LOG=$GIT_STATUSDIR/FDS_log

  MESSAGE_FILE=$GIT_STATUSDIR/message

  MKDIR $GIT_STATUSDIR
# remove untracked files, revert repo files, update to latest revision
  cd $fdsrepo

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
    git remote update
    git merge origin/development
  fi

# get info for smokeview
  cd $SMV_SOURCE
  if [ ! -e $GIT_QT_FILE ]; then
    touch $GIT_QT_FILE
  fi
  THIS_QT_REVISION=`git log --abbrev-commit $QUICKTRIGGER | head -1 | awk '{print $2}'`
  LAST_QT_REVISION=`cat $GIT_QT_FILE`

  THIS_T_REVISION=`git log --abbrev-commit $TRIGGER | head -1 | awk '{print $2}'`
  LAST_T_REVISION=`cat $GIT_T_FILE`

  THIS_SMVAUTHOR=`git log . | head -2 | tail -1 | awk '{print $2}'`
  if [ ! -e $GIT_SMV_FILE ]; then
    touch $GIT_SMV_FILE
  fi
  
  THIS_SMV_REVISION=`git log --abbrev-commit . | head -1 | awk '{print $2}'`
  LAST_SMV_REVISION=`cat $GIT_SMV_FILE`
  git log . | head -5 | tail -1 > $GIT_SMV_LOG

# get info for FDS
  cd $FDS_SOURCE
  THIS_FDSAUTHOR=`git log . | head -2 | tail -1 | awk '{print $2}'`
  if [ ! -e $GIT_FDS_FILE ]; then
    touch $GIT_FDS_FILE
  fi
  THIS_FDS_REVISION=`git log --abbrev-commit . | head -1 | awk '{printf $2}'`
  LAST_FDS_REVISION=`cat $GIT_FDS_FILE`
  git log . | head -5 | tail -1 > $GIT_FDS_LOG

  if [ "$option" == "" ]; then
    if [[ $THIS_SMV_REVISION == $LAST_SMV_REVISION && $THIS_FDS_REVISION == $LAST_FDS_REVISION ]] ; then
      exit
    fi
  else
    if [[ $THIS_QT_REVISION == $LAST_QT_REVISION && $THIS_T_REVISION == $LAST_T_REVISION ]] ; then
      exit
    fi
    if [[ $THIS_QT_REVISION != $LAST_QT_REVISION ]] ; then
      SMOKEBOT_LITE=1
    fi
    if [[ $THIS_T_REVISION != $LAST_T_REVISION ]] ; then
      SMOKEBOT_LITE=
    fi
  fi

  rm -f $MESSAGE_FILE
  if [ "$option" == "" ]; then
    if [[ $THIS_SMV_REVISION != $LAST_SMV_REVISION ]] ; then
      echo $THIS_SMV_REVISION>$GIT_SMV_FILE
      echo -e "smokeview source has changed. $LAST_SMV_REVISION->$THIS_SMV_REVISION($THIS_SMVAUTHOR)" >> $MESSAGE_FILE
      cat $GIT_SMV_LOG >> $MESSAGE_FILE
    fi
    if [[ $THIS_FDS_REVISION != $LAST_FDS_REVISION ]] ; then
      echo $THIS_FDS_REVISION>$GIT_FDS_FILE
      echo -e "FDS source has changed. $LAST_FDS_REVISION->$THIS_FDS_REVISION($THIS_FDSAUTHOR)" >> $MESSAGE_FILE
      cat $GIT_FDS_LOG >> $MESSAGE_FILE
    fi
  else
    if [[ $THIS_QT_REVISION != $LAST_QT_REVISION ]] ; then
      echo $THIS_QT_REVISION>$GIT_QT_FILE
      echo -e "quick trigger file has changed. " >> $MESSAGE_FILE
      cat $GIT_SMV_LOG >> $MESSAGE_FILE
    fi
    if [[ $THIS_T_REVISION != $LAST_T_REVISION ]] ; then
      echo $THIS_T_REVISION>$GIT_T_FILE
      echo -e "trigger file has changed. " >> $MESSAGE_FILE
      cat $GIT_SMV_LOG >> $MESSAGE_FILE
    fi
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
    mkdir -p $DIR
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
   cd $fdsrepo
   chmod -R go+r *
}

clean_repo()
{
  curdir=`pwd`
  dir=$1
  cd $dir
  git clean -dxf &> /dev/null
  git add . &> /dev/null
  git reset --hard HEAD &> /dev/null
  cd $curdir
}


clean_smokebot_history()
{
   
   # Clean Smokebot metafiles
   MKDIR $SMOKEBOT_RUNDIR > /dev/null
   cd $SMOKEBOT_RUNDIR
   MKDIR guides > /dev/null
   MKDIR $HISTORY_DIR > /dev/null
   MKDIR $OUTPUT_DIR > /dev/null
   rm -rf $OUTPUT_DIR/* > /dev/null
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

update_cfast()
{
   cd $SMOKEBOT_HOME_DIR

   # Check to see if CFAST repository exists
   updateclean=
   echo "cfast repo"
   # If yes, then update the CFAST repository and compile CFAST
   if [ -e "$cfastrepo" ] ; then
      cd $cfastrepo
      IS_DIRTY=`git describe --long --dirty | grep dirty | wc -l`
      if [ "$CLEANREPO" == "1" ]; then
        echo "   cleaning"
        if [ "$IS_DIRTY" == "1" ]; then
          echo "The repo $cfastrepo has uncommitted changes"
          echo "Commit or revert these changes or re-run"
          echo "smokebot without the -c (clean) option"
          exit
        fi
        clean_repo $cfastrepo
        updateclean="1"
      fi

      # Update to latest GIT revision
      if [ "$UPDATEREPO" == "1" ]; then
        echo "   updating"
        if [ "$IS_DIRTY" == "1" ]; then
          echo "The repo $cfastrepo has uncommitted changes."
          echo "Commit or revert these changes or re-run"
          echo "smokebot without the -u (update) option"
          exit
        fi
        echo "Updating cfast repo:" >> $OUTPUT_DIR/stage0a
        git remote update >> $OUTPUT_DIR/stage0a 2>&1
        git merge origin/master >> $OUTPUT_DIR/stage0a 2>&1
        updateclean="1"
      fi
      if [ "$updateclean" == "" ]; then
         echo "   not cleaned or updated"
      fi 
   else
      echo "The cfast repo $cfastrepo does not exist"
      echo "Aborting  smokebot"
      exit
   fi

}

compile_cfast()
{
   cd $SMOKEBOT_HOME_DIR

    # Build CFAST
    echo "Building"
    echo "   cfast release"
    cd $cfastrepo/Build/CFAST/${COMPILER}_${platform}${size}
    rm -f cfast7_${platform}${size}
    make --makefile ../makefile clean &> /dev/null
    ./make_cfast.sh >> $OUTPUT_DIR/stage1a 2>&1

   # Check for errors in CFAST compilation
   cd $cfastrepo/Build/CFAST/${COMPILER}_${platform}${size}
   if [ -e "cfast7_${platform}${size}" ]
   then
      stage0_success=true
   else
      echo "Errors from Stage 0 - CFAST:" >> $ERROR_LOG
      echo "CFAST failed to compile" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage1a >> $ERROR_LOG
      echo "" >> $ERROR_LOG
      THIS_CFAST_FAILED=1
   fi
}

#  ============================
#  = Stage 1 - GIT operations =
#  ============================

clean_FDS_repo()
{
   # Check to see if FDS repository exists
   updateclean=
   if [ -e "$fdsrepo" ]
   then
      echo FDS-SMV repo
      if [ "$CLEANREPO" == "1" ]; then
        cd $fdsrepo
        echo "   cleaning"
        IS_DIRTY=`git describe --long --dirty | grep dirty | wc -l`
        if [ "$IS_DIRTY" == "1" ]; then
          echo "The repo $fdsrepo has uncommitted changes."
          echo "Commit or revert these changes or re-run"
          echo "smokebot without the -c (clean) option"
          exit
        fi
        clean_repo $fdsrepo/Verification
        clean_repo $fdsrepo/SMV
        clean_repo $fdsrepo/FDS_Source
        clean_repo $fdsrepo/FDS_Compilation
        clean_repo $fdsrepo/Manuals
        updateclean="1"
      fi
   else
      echo "The FDS repository $fdsrepo does not exist." >> $OUTPUT_DIR/stage0b 2>&1
      echo "Aborting smokebot" >> $OUTPUT_DIR/stage0b 2>&1
      exit
   fi
}

do_FDS_checkout()
{
   cd $fdsrepo
 
   CURRENT_BRANCH=`git rev-parse --abbrev-ref HEAD`
   if [[ "$BRANCH" != "" ]] ; then
     if [[ `git branch | grep $BRANCH` == "" ]] ; then 
        echo "Error: the branch $BRANCH does not exist."
        echo "Aborting smokebot"
        exit
     fi
     if [[ "$BRANCH" != "$CURRENT_BRANCH" ]] ; then
        echo "Checking out branch $BRANCH." >> $OUTPUT_DIR/stage0b 2>&1
        git checkout $BRANCH
     fi
   else
      BRANCH=$CURRENT_BRANCH
   fi
   if [ "$UPDATEREPO" == "1" ]; then
     echo "   updating"
     IS_DIRTY=`git describe --long --dirty | grep dirty | wc -l`
     if [ "$IS_DIRTY" == "1" ]; then
       echo "The repo $fdsrepo has uncommitted changes."
       echo "Commit or revert these changes or re-run"
       echo "smokebot without the -u (update) option"
       exit
     fi
     echo "Updating branch $BRANCH." >> $OUTPUT_DIR/stage0b 2>&1
     git remote update >> $OUTPUT_DIR/stage0b 2>&1
     git merge origin/development >> $OUTPUT_DIR/stage0b 2>&1
     echo "Updating submodules." >> $OUTPUT_DIR/stage0b 2>&1
     git submodule foreach git remote update >> $OUTPUT_DIR/stage0b 2>&1
     git submodule foreach git merge origin/master  >> $OUTPUT_DIR/stage0b 2>&1
     updateclean="1"
   fi
   if [ "$updateclean" == "" ]; then
      echo "   not cleaned or updated"
   fi 
   GIT_REVISION=`git describe --long --dirty`
   GIT_SHORTHASH=`git rev-parse --short HEAD`
   GIT_LONGHASH=`git rev-parse HEAD`
   GIT_DATE=`git log -1 --format=%cd --date=local $GIT_SHORTHASH`
}

check_FDS_checkout()
{
   cd $fdsrepo
   # Check for GIT errors
   stage0b_success=true
}

#  ==================================
#  = Stage 1a/b - Compile FDS debug =
#  ==================================

compile_fds_mpi_db()
{
   # Clean and compile mpi FDS debug
   echo "   FDS"
   echo "      debug"
   cd $fdsrepo/FDS_Compilation/mpi_${COMPILER}_${platform}${size}$IB$DB
   rm -f fds_mpi_${COMPILER}_${platform}${size}$IB$DB
   make --makefile ../makefile clean &> /dev/null
   ./make_fds.sh &> $OUTPUT_DIR/stage1b
}

check_compile_fds_mpi_db()
{
   # Check for errors in FDS debug compilation
   cd $fdsrepo/FDS_Compilation/mpi_${COMPILER}_${platform}${size}$IB$DB
   if [ -e "fds_mpi_${COMPILER}_${platform}${size}$IB$DB" ]
   then
      stage1b_fdsdb_success=true
   else
      echo "Errors from Stage 1b - Compile FDS MPI debug:" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage1b >> $ERROR_LOG
      echo "" >> $ERROR_LOG
      THIS_FDS_FAILED=1
   fi

   # Check for compiler warnings/remarks
   if [[ `grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage1b| grep -v 'find atom' | grep -v 'feupdateenv is not implemented'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 1b warnings:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage1b | grep -v 'find atom' | grep -v 'feupdateenv is not implemented'>> $WARNING_LOG
      echo "" >> $WARNING_LOG
   # if the executable does not exist then an email has already been sent
      if [ -e "fds_mpi_${COMPILER}_${platform}${size}$IB$DB" ] ; then
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
   echo "Running verification cases"
   if [ "$CLEANREPO" == "1" ]; then
     echo "   cleaning"
     cd $fdsrepo/Verification
     clean_repo $fdsrepo/Verification
   fi

   #  =====================
   #  = Run all SMV cases =
   #  =====================

   echo "   debug"
   cd $fdsrepo/Verification/scripts

   # Submit SMV verification cases and wait for them to start
   echo 'Running SMV verification cases:' >> $OUTPUT_DIR/stage3a 2>&1
   ./Run_SMV_Cases.sh -c $cfastrepo -I $COMPILER $USEINSTALL2 -m 2 -d -q $SMOKEBOT_QUEUE -j $JOBPREFIX >> $OUTPUT_DIR/stage3a 2>&1

   # Wait for SMV verification cases to end
   wait_verification_cases_debug_end

}

check_verification_cases_debug()
{
   # Scan and report any errors in FDS verification cases
   cd $fdsrepo/Verification

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
#  = Stage 1c - Compile FDS release =
#  ====================================

compile_fds_mpi()
{
   # Clean and compile FDS
   echo "      release"
   cd $fdsrepo/FDS_Compilation/mpi_${COMPILER}_${platform}${size}$IB
   rm -f fds_mpi_${COMPILER}_${platform}${size}$IB
   make --makefile ../makefile clean &> /dev/null
   ./make_fds.sh &> $OUTPUT_DIR/stage1c
}

check_compile_fds_mpi()
{
   # Check for errors in FDS compilation
   cd $fdsrepo/FDS_Compilation/mpi_${COMPILER}_${platform}${size}$IB
   if [ -e "fds_mpi_${COMPILER}_${platform}${size}$IB" ]
   then
      stage1c_fdsrel_success=true
   else
      echo "Errors from Stage 1c - Compile FDS release:" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage1c >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   # 'performing multi-file optimizations' and 'generating object file' are part of a normal compile
   if [[ `grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage1c | grep -v 'find atom' | grep -v 'performing multi-file optimizations' | grep -v 'generating object file'| grep -v 'feupdateenv is not implemented'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 1c warnings:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage1c | grep -v 'find atom' | grep -v 'performing multi-file optimizations' | grep -v 'generating object file'| grep -v 'feupdateenv is not implemented' >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

#  ======================================
#  = Stage 2a - Compile SMV utilities =
#  ======================================

compile_smv_utilities()
{
   echo "   smokeview utilities"
   echo "" > $OUTPUT_DIR/stage2a
   if [ "$haveCC" == "1" ] ; then
   if [ "$SSH" == "" ] ; then 
   # smokeview libraries
   echo "      libraries"
   cd $fdsrepo/SMV/Build/LIBS/lib_${platform}_${COMPILER}${size}
   echo 'Building Smokeview libraries:' >> $OUTPUT_DIR/stage2a 2>&1
   ./makelibs.sh >> $OUTPUT_DIR/stage2a 2>&1

   # smokezip:
   echo "      smokezip"
   cd $fdsrepo/SMV/Build/smokezip/${COMPILER}_${platform}${size}
   rm -f *.o smokezip_${platform}${size}
   echo 'Compiling smokezip:' >> $OUTPUT_DIR/stage2a 2>&1
   ./make_smokezip.sh >> $OUTPUT_DIR/stage2a 2>&1
   echo "" >> $OUTPUT_DIR/stage2a 2>&1
   
   # smokediff:
   echo "      smokediff"
   cd $fdsrepo/SMV/Build/smokediff/${COMPILER}_${platform}${size}
   rm -f *.o smokediff_${platform}${size}
   echo 'Compiling smokediff:' >> $OUTPUT_DIR/stage2a 2>&1
   ./make_smokediff.sh >> $OUTPUT_DIR/stage2a 2>&1
   echo "" >> $OUTPUT_DIR/stage2a 2>&1
   
   # background
   echo "      background"
   cd $fdsrepo/SMV/Build/background/${COMPILER}_${platform}${size}
   rm -f *.o background
   echo 'Compiling background:' >> $OUTPUT_DIR/stage2a 2>&1
   ./make_background.sh >> $OUTPUT_DIR/stage2a 2>&1
   
  # wind2fds:
   echo "      wind2fds"
   cd $fdsrepo/SMV/Build/wind2fds/${COMPILER}_${platform}${size}
   rm -f *.o wind2fds_${platform}${size}
   echo 'Compiling wind2fds:' >> $OUTPUT_DIR/stage2a 2>&1
   ./make_wind2fds.sh >> $OUTPUT_DIR/stage2a 2>&1
   echo "" >> $OUTPUT_DIR/stage2a 2>&1
   else
   $SSH \( \
   cd $fdsrepo/SMV/Build/LIBS/lib_${platform}_${COMPILER}${size} \; \
   echo 'Building Smokeview libraries:' >> $OUTPUT_DIR/stage2a 2>&1 \; \
   ./makelibs.sh >> $OUTPUT_DIR/stage2a 2>&1 \; \
   cd $fdsrepo/SMV/Build/smokezip/${COMPILER}_${platform}${size} \; \
   rm -f *.o smokezip_${platform}${size} \; \
   echo 'Compiling smokezip:' >> $OUTPUT_DIR/stage2a 2>&1 \; \
   ./make_smokezip.sh >> $OUTPUT_DIR/stage2a 2>&1 \; \
   echo "" >> $OUTPUT_DIR/stage2a 2>&1 \; \
   cd $fdsrepo/SMV/Build/smokediff/${COMPILER}_${platform}${size} \; \
   rm -f *.o smokediff_${platform}${size} \; \
   echo 'Compiling smokediff:' >> $OUTPUT_DIR/stage2a 2>&1 \; \
   ./make_smokediff.sh >> $OUTPUT_DIR/stage2a 2>&1 \; \
   echo "" >> $OUTPUT_DIR/stage2a 2>&1 \; \
   cd $fdsrepo/SMV/Build/background/${COMPILER}_${platform}${size} \; \
   rm -f *.o background \; \
   echo 'Compiling background:' >> $OUTPUT_DIR/stage2a 2>&1 \; \
   ./make_background.sh >> $OUTPUT_DIR/stage2a 2>&1 \; \
   cd $fdsrepo/SMV/Build/wind2fds/${COMPILER}_${platform}${size} \; \
   echo 'Compiling wind2fds:' >> $OUTPUT_DIR/stage2a 2>&1 \; \
   ./make_wind2fds.sh >> $OUTPUT_DIR/stage2a 2>&1 \; \
   echo "" >> $OUTPUT_DIR/stage2a 2>&1  \)
   fi
   else
   echo "Warning: smokeview and utilities not built - C compiler not available" >> $OUTPUT_DIR/stage2a 2>&1
   fi
}

is_file_installed()
{
  program=$1
  notfound=`$program -help | tail -1 | grep "not found" | wc -l`
  if [ "$notfound" == "1" ] ; then
    stage2a_success="0"
    echo "***error: $program not installed" >> $OUTPUT_DIR/stage2a
  fi
}

check_smv_utilities()
{
   if [ "$haveCC" == "1" ] ; then
     # Check for errors in SMV utilities compilation
     cd $fdsrepo
     if [ -e "$fdsrepo/SMV/Build/smokezip/${COMPILER}_${platform}${size}/smokezip_${platform}${size}" ]  && \
        [ -e "$fdsrepo/SMV/Build/smokediff/${COMPILER}_${platform}${size}/smokediff_${platform}${size}" ]  && \
        [ -e "$fdsrepo/SMV/Build/wind2fds/${COMPILER}_${platform}${size}/wind2fds_${platform}${size}" ]  && \
        [ -e "$fdsrepo/SMV/Build/background/${COMPILER}_${platform}${size}/background" ]
     then
        stage2a_success="1"
     else
        stage2a_success="0"
        echo "Errors from Stage 2a - Compile SMV utilities:" >> $ERROR_LOG
        cat $OUTPUT_DIR/stage2a >> $ERROR_LOG
        echo "" >> $ERROR_LOG
     fi
   else
     stage2a_success="1"
     is_file_installed smokeview
     is_file_installed smokezip
     is_file_installed smokediff
     is_file_installed wind2fds
     is_file_installed background
     if [ "$stage2a_success" == "0" ] ; then
        echo "Errors from Stage 2a - Smokeview and utilities:" >> $ERROR_LOG
        stage2a_success="1"
        cat $OUTPUT_DIR/stage2a >> $ERROR_LOG
        echo "" >> $ERROR_LOG
     fi
   fi
}

#  ===================================================
#  = Stage 3 - Run verification cases (release mode) =
#  ===================================================

wait_verification_cases_release_end()
{
   # Scans qstat and waits for verification cases to end
   if [[ "$SMOKEBOT_QUEUE" == "none" ]]
   then
     while [[ `ps -u $USER -f | fgrep .fds | grep -v grep` != '' ]]; do
        JOBS_REMAINING=`ps -u $USER -f | fgrep .fds | grep -v grep | wc -l`
        echo "Waiting for ${JOBS_REMAINING} verification cases to complete." >> $OUTPUT_DIR/stage3b
        TIME_LIMIT_STAGE="5"
        check_time_limit
        sleep 60
     done
   else
     while [[ `qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX` != '' ]]; do
        JOBS_REMAINING=`qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX | wc -l`
        echo "Waiting for ${JOBS_REMAINING} verification cases to complete." >> $OUTPUT_DIR/stage3b
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
     echo "   clean"
     cd $fdsrepo/Verification
     clean_repo $fdsrepo/Verification
   fi
   echo "   release"
   # Start running all SMV verification cases
   cd $fdsrepo/Verification/scripts
   echo 'Running SMV verification cases:' >> $OUTPUT_DIR/stage3b 2>&1
   ./Run_SMV_Cases.sh -c $cfastrepo -I $COMPILER $USEINSTALL2 $RUN_OPENMP -q $SMOKEBOT_QUEUE -j $JOBPREFIX >> $OUTPUT_DIR/stage3b 2>&1

   # Wait for all verification cases to end
   wait_verification_cases_release_end
}

check_verification_cases_release()
{
   # Scan and report any errors in FDS verification cases
   cd $fdsrepo/Verification

   if [[ `grep -rIi 'Run aborted' $OUTPUT_DIR/stage3b` == "" ]] && \
      [[ `grep -rIi 'Segmentation' Visualization/* WUI/* Immersed_Boundary_Method/* ` == "" ]] && \
      [[ `grep -rI 'ERROR:' Visualization/* WUI/* Immersed_Boundary_Method/* ` == "" ]] && \
      [[ `grep -rIi 'STOP: Numerical' Visualization/* WUI/* Immersed_Boundary_Method/* ` == "" ]] && \
      [[ `grep -rIi  -A 20 'forrtl' Visualization/* WUI/* Immersed_Boundary_Method/* ` == "" ]]
   then
      stage3b_success=true
   else
      grep -rIi 'Run aborted' $OUTPUT_DIR/stage3b > $OUTPUT_DIR/stage3b_errors
      grep -rIi 'Segmentation' Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage3b_errors
      grep -rI 'ERROR:' Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage3b_errors
      grep -rIi 'STOP: Numerical' Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage3b_errors
      grep -rIi -A 20 'forrtl' Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage3b_errors

      echo "Errors from Stage 3b - Run verification cases (release mode):" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage3b_errors >> $ERROR_LOG
      echo "" >> $ERROR_LOG
      THIS_FDS_FAILED=1
   fi

      
   if [[ `grep 'Warning' -rI $OUTPUT_DIR/stage3b` == "" ]] 
   then
      no_warnings=true
   else
      echo "Stage 3b warnings:" >> $WARNING_LOG
      grep 'Warning' -rI $OUTPUT_DIR/stage3b >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

#  ================================
#  = Stage 2b - Compile SMV debug =
#  ================================

compile_smv_db()
{
   if [ "$haveCC" == "1" ] ; then
   if [ "$SSH" == "" ] ; then
   # Clean and compile SMV debug
   echo "   smokeview"
   echo "      debug"
   cd $fdsrepo/SMV/Build/smokeview/${COMPILER}_${platform}${size}
   rm -f smokeview_${platform}${size}_db
   ./make_smv_db.sh &> $OUTPUT_DIR/stage2b
   else
   $SSH \(
   cd $fdsrepo/SMV/Build/smokeview/${COMPILER}_${platform}${size} \; \
   rm -f smokeview_${platform}${size}_db \; \
   ./make_smv_db.sh &> $OUTPUT_DIR/stage2b \)
   fi
   fi
}

check_compile_smv_db()
{
   if [ "$haveCC" == "1" ] ; then
   # Check for errors in SMV debug compilation
   cd $fdsrepo/SMV/Build/smokeview/${COMPILER}_${platform}${size}
   if [ -e "smokeview_${platform}${size}_db" ]
   then
      stage2b_success=true
   else
      echo "Errors from Stage 2b - Compile SMV debug:" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage2b >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   # grep -v 'feupdateenv ...' ignores a known FDS MPI compiler warning (http://software.intel.com/en-us/forums/showthread.php?t=62806)
   if [[ `grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage2b | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 6a warnings:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage2b | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked' >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
   fi
}

#  =============================================
#  = Stage 4a - Make SMV pictures (debug mode) =
#  =============================================

make_smv_pictures_db()
{
   # Run Make SMV Pictures script (debug mode)
   if [ "$SSH" == "" ]; then
   echo "making smokeview images"
   cd $fdsrepo/Verification/scripts
   ./Make_SMV_Pictures.sh $USEINSTALL -d 2>&1 &> $OUTPUT_DIR/stage4a_orig
   grep -v FreeFontPath $OUTPUT_DIR/stage4a_orig > $OUTPUT_DIR/stage4a
   else
   $SSH \( cd $fdsrepo/Verification/scripts \; \
   ./Make_SMV_Pictures.sh $USEINSTALL -d 2>&1 \| grep -v FreeFontPath &> $OUTPUT_DIR/stage4a \)
   fi
}

check_smv_pictures_db()
{
   # Scan and report any errors in make SMV pictures process
   echo "   checking"
   cd $SMOKEBOT_RUNDIR
   if [[ `grep -I -E "Segmentation|Error" $OUTPUT_DIR/stage4a` == "" ]]
   then
      stage4a_success=true
   else
      cp $OUTPUT_DIR/stage4a $OUTPUT_DIR/stage4a_errors

      echo "Errors from Stage 4a - Make SMV pictures (debug mode):" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage4a_errors >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Scan for and report any warnings in make SMV pictures process
   cd $SMOKEBOT_RUNDIR
   if [[ `grep -I -E "Warning" $OUTPUT_DIR/stage4a` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 6b - Make SMV pictures (debug mode):" >> $WARNING_LOG
      grep -I -E "Warning" $OUTPUT_DIR/stage4a >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

#  ==================================
#  = Stage 2c - Compile SMV release =
#  ==================================

compile_smv()
{
   if [ "$haveCC" == "1" ] ; then
   if [ "$SSH" == "" ] ; then
   # Clean and compile SMV
   echo "      release"
   cd $fdsrepo/SMV/Build/smokeview/${COMPILER}_${platform}${size}
   rm -f smokeview_${platform}${size}
   ./make_smv.sh $TESTFLAG &> $OUTPUT_DIR/stage2c
   else
   $SSH \( \
   cd $fdsrepo/SMV/Build/smokeview/${COMPILER}_${platform}${size} \; \
   rm -f smokeview_${platform}${size} \; \
   ./make_smv.sh $TESTFLAG &> $OUTPUT_DIR/stage2c \)
   fi
   fi
}

check_compile_smv()
{
   if [ "$haveCC" == "1" ] ; then
   # Check for errors in SMV release compilation
   cd $fdsrepo/SMV/Build/smokeview/${COMPILER}_${platform}${size}
   if [ -e "smokeview_${platform}${size}" ]
   then
      stage2c_smv_success=true
   else
      echo "Errors from Stage 2c - Compile SMV release:" >> $ERROR_LOG
      echo "The program smokeview_${platform}${size} does not exist."
      cat $OUTPUT_DIR/stage2c >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   # grep -v 'feupdateenv ...' ignores a known FDS MPI compiler warning (http://software.intel.com/en-us/forums/showthread.php?t=62806)
   if [[ `grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage2c | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 6c warnings:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage2c | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked' >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
   fi
}

#  ===============================================
#  = Stage 4 - Make SMV pictures (release mode) =
#  ===============================================

make_smv_pictures()
{
   # Run Make SMV Pictures script (release mode)
   echo Generating images 
   if [ "$SSH" == "" ]; then
   cd $fdsrepo/Verification/scripts
   ./Make_SMV_Pictures.sh $TESTFLAG $USEINSTALL 2>&1 &> $OUTPUT_DIR/stage4b_orig
   grep -v FreeFontPath $OUTPUT_DIR/stage4b_orig &> $OUTPUT_DIR/stage4b
   else
   $SSH \( cd $fdsrepo/Verification/scripts \; \
   ./Make_SMV_Pictures.sh $TESTFLAG $USEINSTALL 2>&1 \| grep -v FreeFontPath &> $OUTPUT_DIR/stage4b \)
   fi
}

check_smv_pictures()
{
   # Scan and report any errors in make SMV pictures process
   cd $SMOKEBOT_RUNDIR
   echo "   checking"
   if [[ `grep -I -E "Segmentation|Error" $OUTPUT_DIR/stage4b` == "" ]]
   then
      stage4b_smvpics_success=true
   else
      cp $OUTPUT_DIR/stage4b  $OUTPUT_DIR/stage4b_errors

      echo "Errors from Stage 4 - Make SMV pictures (release mode):" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage4b >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi
   if [ ! "$web_DIR" == "" ]; then
     if [ -d "$WEBFROMDIR" ]; then
       CURDIR=`pwd`
       cd $web_DIR
       rm -rf images
       rm -rf manuals
       rm index.html
       cd $WEBFROMDIR
       cp -r * $web_DIR/.
       cd $CURDIR
     fi
   fi

}

#  ===============================================
#  = Stage 4c - Make SMV movies (release mode) =
#  ===============================================

make_smv_movies()
{
   cd $fdsrepo/Verification
   scripts/Make_SMV_Movies.sh 2>&1  &> $OUTPUT_DIR/stage4c
}

check_smv_movies()
{
   cd $SMOKEBOT_RUNDIR
   echo make smokeview movies
   if [[ `grep -I -E "Segmentation|Error" $OUTPUT_DIR/stage4c` == "" ]]
   then
      stage4c_success=true
   else
      cp $OUTPUT_DIR/stage4c  $OUTPUT_DIR/stage4c_errors

      echo "Errors from Stage 4c - Make SMV movies " >> $ERROR_LOG
      cat $OUTPUT_DIR/stage4c >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Scan for and report any warnings in make SMV pictures process
   cd $SMOKEBOT_RUNDIR
   if [[ `grep -I -E "Warning" $OUTPUT_DIR/stage4c` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 6e - Make SMV movies (release mode):" >> $WARNING_LOG
      grep -I -E "Warning" $OUTPUT_DIR/stage4c >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
   if [ ! "$web_DIR" == "" ]; then
     if [ -d "$WEBFROMDIR" ]; then 
       CURDIR=`pwd`
       cd $web_DIR
       rm -rf *
       cd $WEBFROMDIR
       cp -r * $web_DIR/.
       cd $CURDIR
     fi
   fi

}

#  ======================================
#  = FDS run time statistics =
#  ======================================

generate_timing_stats()
{
   echo "Timing stats"
   echo "   generating"
   cd $fdsrepo/Verification/scripts/
   export QFDS="$fdsrepo/Verification/scripts/copyout.sh"
   export RUNCFAST="$fdsrepo/Verification/scripts/copyout.sh"
   export RUNTFDS="$fdsrepo/Verification/scripts/copyout.sh"

   cd $fdsrepo/Verification
   scripts/SMV_Cases.sh
   scripts/GEOM_Cases.sh

   cd $fdsrepo/Utilities/Scripts
   ./fds_timing_stats.sh smokebot > smv_timing_stats.csv
   cd $fdsrepo/Utilities/Scripts
   ./fds_timing_stats.sh smokebot 1 > smv_benchmarktiming_stats.csv
   TOTAL_SMV_TIMES=`tail -1 smv_benchmarktiming_stats.csv`
}

archive_timing_stats()
{
  echo "   archiving"
  cd $fdsrepo/Utilities/Scripts
  cp smv_timing_stats.csv "$HISTORY_DIR/${GIT_REVISION}_timing.csv"
  cp smv_benchmarktiming_stats.csv "$HISTORY_DIR/${GIT_REVISION}_benchmarktiming.csv"
  TOTAL_SMV_TIMES=`tail -1 smv_benchmarktiming_stats.csv`
  if [ "$UPLOADRESULTS" == "1" ]; then
    cd $fdsrepo/Utilities/Firebot
    ./smvstatus_updatepub.sh -F
  fi
  if [ ! "$web_DIR" == "" ]; then
    cd $fdsrepo/Utilities/Firebot
    ./make_smv_summary.sh > $web_DIR/index.html
  fi
}

#  ===================================
#  = Stage 5 - Build smokview guides =
#  ===================================

check_guide()
{
   stage=$1
   directory=$2
   document=$3
   label=$4

   # Scan and report any errors in build process for guides

   SMOKEBOT_MANDIR=
   if [ ! "$web_DIR" == "" ]; then
     if [ -d $web_DIR/manuals ]; then
       SMOKEBOT_MANDIR=$web_DIR/manuals
     fi
   fi

   cd $SMOKEBOT_RUNDIR
   if [[ `grep "! LaTeX Error:" -I $stage` == "" ]]; then
     if [ ! "$SMOKEBOT_MANDIR" == "" ]; then
       cp $directory/$document $SMOKEBOT_MANDIR/.
     fi
     if [ -d $SMV_SUMMARY/manuals ] ; then
       cp $directory/$document $SMV_SUMMARY/manuals/.
     fi
     cp $directory/$document $NEWGUIDE_DIR/.
     chmod 664 $NEWGUIDE_DIR/$document
   else
      echo "Errors from Stage 5 - Build FDS-SMV Guides:" >> $ERROR_LOG
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
      echo "Stage 5 warnings:" >> $WARNING_LOG
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
  
   ./make_guide.sh &> $OUTPUT_DIR/stage5_$document

   # Check guide for completion and copy to website if successful
   check_guide $OUTPUT_DIR/stage5_$document $directory $document.pdf $label
}

#  =====================================================
#  = Build status reporting - email and save functions =
#  =====================================================

save_build_status()
{
   STOP_TIME=$(date)
   STOP_TIME_INT=$(date +%s)
   cd $SMOKEBOT_RUNDIR
   # Save status outcome of build to a text file
   if [[ -e $WARNING_LOG && -e $ERROR_LOG ]]
   then
     echo "***Warnings:" >> $ERROR_LOG
     cat $WARNING_LOG >> $ERROR_LOG
     echo "   build failure and warnings for version: ${GIT_REVISION}, branch: $BRANCH."
     echo "Build failure and warnings;$GIT_DATE;$GIT_SHORTHASH;$GIT_LONGHASH;${GIT_REVISION};$BRANCH;$STOP_TIME_INT;3;$TOTAL_SMV_TIMES" > "$HISTORY_DIR/${GIT_REVISION}.txt"
     cat $ERROR_LOG > "$HISTORY_DIR/${GIT_REVISION}_errors.txt"

   # Check for errors only
   elif [ -e $ERROR_LOG ]
   then
      echo "   build failure for version: ${GIT_REVISION}, branch: $BRANCH."
      echo "Build failure;$GIT_DATE;$GIT_SHORTHASH;$GIT_LONGHASH;${GIT_REVISION};$BRANCH;$STOP_TIME_INT;3;$TOTAL_SMV_TIMES" > "$HISTORY_DIR/${GIT_REVISION}.txt"
      cat $ERROR_LOG > "$HISTORY_DIR/${GIT_REVISION}_errors.txt"

   # Check for warnings only
   elif [ -e $WARNING_LOG ]
   then
      echo "   build success with warnings for version: ${GIT_REVISION}, branch: $BRANCH."
      echo "Build success with warnings;$GIT_DATE;$GIT_SHORTHASH;$GIT_LONGHASH;${GIT_REVISION};$BRANCH;$STOP_TIME_INT;2;$TOTAL_SMV_TIMES" > "$HISTORY_DIR/${GIT_REVISION}.txt"
      cat $WARNING_LOG > "$HISTORY_DIR/${GIT_REVISION}_warnings.txt"

   # No errors or warnings
   else
      echo "   build success for version: ${GIT_REVISION}, branch: $BRANCH."
      echo "Build success!;$GIT_DATE;$GIT_SHORTHASH;$GIT_LONGHASH;${GIT_REVISION};$BRANCH;$STOP_TIME_INT;1;$TOTAL_SMV_TIMES" > "$HISTORY_DIR/${GIT_REVISION}.txt"
   fi
}

email_build_status()
{
   if [[ "$THIS_FDS_FAILED" == "1" ]] ; then
     mailTo="$mailToFDS"
   fi
   if [[ "$THIS_CFAST_FAILED" == "1" ]] ; then
     mailTo="$mailToCFAST"
   fi
   if [[ "$MAILTO" != "" ]] ; then
     mailTo="$MAILTO"
   fi
   echo $THIS_FDS_FAILED>$FDS_STATUS_FILE
   stop_time=`date`
   echo "----------------------------------------------" > $TIME_LOG
   echo "             host: $hostname " >> $TIME_LOG
   echo "            start: $start_time " >> $TIME_LOG
   echo "             stop: $stop_time " >> $TIME_LOG
   echo "            setup: $DIFF_PRELIM" >> $TIME_LOG
   echo "   build software: $DIFF_BUILDSOFTWARE" >> $TIME_LOG
   echo "        run cases: $DIFF_RUNCASES" >> $TIME_LOG
   echo "    make pictures: $DIFF_MAKEPICTURES" >> $TIME_LOG
if [ "$MAKEMOVIES" == "1" ]; then
   echo "      make movies: $DIFF_MAKEMOVIES" >> $TIME_LOG
fi
   echo "      make guides: $DIFF_MAKEGUIDES" >> $TIME_LOG
   echo "            total: $DIFF_SCRIPT_TIME" >> $TIME_LOG
   echo "benchmark time(s): $TOTAL_SMV_TIMES" >> $TIME_LOG
if [ "$RUNAUTO" == "y" ]; then
   echo "FDS revisions: old: $LAST_FDS_REVISION new: $THIS_FDS_REVISION" >> $TIME_LOG
   echo "SMV revisions: old: $LAST_SMV_REVISION new: $THIS_SMV_REVISION" >> $TIME_LOG
fi
if [ "$RUNAUTO" == "Y" ]; then
   echo "FDS revisions: $THIS_SMV_REVISION" >> $TIME_LOG
   echo "SMV revisions: $THIS_FDS_REVISION" >> $TIME_LOG
fi
if [ "$RUNAUTO" == "" ]; then
   echo "SMV revisions: $THIS_SMV_REVISION" >> $TIME_LOG
fi
  if [[ $THIS_SMV_REVISION != $LAST_SMV_REVISION ]] ; then
    cat $GIT_SMV_LOG >> $TIME_LOG
  fi
  if [[ $THIS_FDS_REVISION != $LAST_FDS_REVISION ]] ; then
    cat $GIT_FDS_LOG >> $TIME_LOG
  fi
   cd $SMOKEBOT_RUNDIR
   # Check for warnings and errors
   if [ ! "$WEB_URL" == "" ]; then
     echo "Smokebot summary: $WEB_URL" >> $TIME_LOG
   fi
   if [ "$UPLOADRESULTS" == "1" ]; then
     echo " Smokebot status: https://goo.gl/gKVSDZ" >> $TIME_LOG
   fi
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
        $UploadGuides $NEWGUIDE_DIR $fdsrepo/Manuals &> /dev/null
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
if [[ $RUNAUTO == "Y" ]] ; then
  run_auto trigger
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

### Stage 0 repo operatoins ###
update_cfast

clean_FDS_repo
do_FDS_checkout
check_FDS_checkout
PRELIM_end=`GET_TIME`
DIFF_PRELIM=`GET_DURATION $PRELIM_beg $PRELIM_end`
echo "Preliminary: $DIFF_PRELIM" >> $STAGE_STATUS

### Stage 1 build cfast and FDS ###
BUILDSOFTWARE_beg=`GET_TIME`
compile_cfast
compile_fds_mpi_db
check_compile_fds_mpi_db

if [ "$SMOKEBOT_LITE" == "" ]; then
if [[ $stage1b_fdsdb_success ]] ; then
   compile_fds_mpi
   check_compile_fds_mpi
fi
fi

### Stage 2 build smokeview ###
compile_smv_utilities
check_smv_utilities

compile_smv_db
check_compile_smv_db

if [ "$SMOKEBOT_LITE" == "" ]; then
  compile_smv
  check_compile_smv
fi

BUILDSOFTWARE_end=`GET_TIME`
DIFF_BUILDSOFTWARE=`GET_DURATION $BUILDSOFTWARE_beg $BUILDSOFTWARE_end`
echo "Build Software: $DIFF_BUILDSOFTWARE" >> $STAGE_STATUS

### Stage 3 run verification cases ###
RUNCASES_beg=`GET_TIME`
if [[ $stage1b_fdsdb_success && "$RUNDEBUG" == "1" ]] ; then
   run_verification_cases_debug
   check_verification_cases_debug
fi

if [ "$SMOKEBOT_LITE" == "" ]; then
  if [[ $stage1c_fdsrel_success ]] ; then
     run_verification_cases_release
     check_verification_cases_release
  fi
fi
RUNCASES_end=`GET_TIME`
DIFF_RUNCASES=`GET_DURATION $RUNCASES_beg $RUNCASES_end`
echo "Run cases: $DIFF_RUNCASES" >> $STAGE_STATUS

### Stage 4 generate images ###
MAKEPICTURES_beg=`GET_TIME`
if [ "$SMOKEBOT_LITE" == "" ]; then
  if [[ $stage1c_fdsrel_success && $stage2c_smv_success ]] ; then
    make_smv_pictures
    check_smv_pictures
  fi
fi
MAKEPICTURES_end=`GET_TIME`
DIFF_MAKEPICTURES=`GET_DURATION $MAKEPICTURES_beg $MAKEPICTURES_end`
echo "Make pictures: $DIFF_MAKEPICTURES" >> $STAGE_STATUS

if [ "$SMOKEBOT_LITE" == "" ]; then
  if [ "$MAKEMOVIES" == "1" ]; then
    MAKEMOVIES_beg=`GET_TIME`
 
    make_smv_movies
    check_smv_movies

    MAKEMOVIES_end=`GET_TIME`
    DIFF_MAKEMOVIES=`GET_DURATION $MAKEMOVIES_beg $MAKEMOVIES_end`
    echo "Make movies: $DIFF_MAKEMOVIES" >> $STAGE_STATUS
fi
fi

if [ "$SMOKEBOT_LITE" == "" ]; then
  if [[ $stage1c_fdsrel_success ]] ; then
    generate_timing_stats
  fi
fi

### Stage 5 build documents ###
MAKEGUIDES_beg=`GET_TIME`
if [ "$SMOKEBOT_LITE" == "" ]; then
  if [[ $stage1c_fdsrel_success && $stage4b_smvpics_success ]] ; then
     echo Making guides
#   echo "   geometry notes"
#  make_guide geom_notes $fdsrepo/Manuals/FDS_User_Guide 'geometry notes'
     echo "   user"
    make_guide SMV_User_Guide $fdsrepo/Manuals/SMV_User_Guide 'SMV User Guide'
     echo "   technical"
    make_guide SMV_Technical_Reference_Guide $fdsrepo/Manuals/SMV_Technical_Reference_Guide 'SMV Technical Reference Guide'
     echo "   verification"
    make_guide SMV_Verification_Guide $fdsrepo/Manuals/SMV_Verification_Guide 'SMV Verification Guide'
  else
    echo Errors found, not building guides
  fi
fi
MAKEGUIDES_end=`GET_TIME`
DIFF_MAKEGUIDES=`GET_DURATION $MAKEGUIDES_beg $MAKEGUIDES_end`
echo "Make guides: $DIFF_MAKEGUIDES" >> $STAGE_STATUS

SCRIPT_TIME_end=`GET_TIME`
DIFF_SCRIPT_TIME=`GET_DURATION $SCRIPT_TIME_beg $SCRIPT_TIME_end`
echo "Total time: $DIFF_SCRIPT_TIME" >> $STAGE_STATUS

### Report results ###
echo Reporting results
set_files_world_readable
save_build_status
 
if [ "$SMOKEBOT_LITE" == "" ]; then
  if [[ $stage1c_fdsrel_success ]] ; then
    archive_timing_stats
  fi
fi
echo "   emailing results"
email_build_status
