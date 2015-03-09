#!/bin/bash

# Smokebot
# This script is a simplified version of Kris Overholt's firebot script.
# It runs the smokeview verification suite (not FDS) on the latest
# revision of the repository.  It does not erase files that are not
# in the repository.  This allows one to test working files before they
# have been committed.  

#  ===================
#  = Input variables =
#  ===================

SMOKEBOT_USERNAME="smokebot"
SMOKEBOT_QUEUE=smokebot
MAKEMOVIES=
RUNAUTO=
BUILDBUNDLE=
RUNDEBUG="1"
OPENMP=
RUN_OPENMP=
TESTFLAG=
RUNFDSCASES=

WEBHOSTNAME=blaze.nist.gov
if [ "$SMOKEBOT_HOSTNAME" != "" ] ; then
WEBHOSTNAME=$SMOKEBOT_HOSTNAME
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

while getopts 'abmo:q:s' OPTION
do
case $OPTION in
  a)
   RUNAUTO="y"
   ;;
  b)
   BUILDBUNDLE="y"
   ;;
  m)
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
  s)
   RUNDEBUG="0"
   ;;
  t)
   TESTFLAG="-t"
   ;;
esac
done
shift $(($OPTIND-1))

if [[ `whoami` == "$SMOKEBOT_USERNAME" ]];
   then
      # Continue along
      :
   else
      echo "Warning: You are running the Smokebot script as an end user."
      echo "This script will definitely modify and/or erase your repository."
      echo "If you wish to continue, run Firebot with the -y option at your own risk."
      echo "Terminating script."
      exit
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

cd
SMOKEBOT_HOME_DIR="`pwd`"
SMOKEBOT_DIR="$SMOKEBOT_HOME_DIR/smokebot"
OUTPUT_DIR="$SMOKEBOT_DIR/output"
export FDS_SVNROOT="$SMOKEBOT_HOME_DIR/FDS-SMV"
export SMV_Summary="$FDS_SVNROOT/Manuals/SMV_Summary"
CFAST_SVNROOT="$SMOKEBOT_HOME_DIR/cfast"
ERROR_LOG=$OUTPUT_DIR/errors
TIME_LOG=$OUTPUT_DIR/timings
WARNING_LOG=$OUTPUT_DIR/warnings
GUIDE_DIR=$SMOKEBOT_DIR/guides
STAGE_STATUS=$OUTPUT_DIR/stage_status
SMV_VG_GUIDE=$FDS_SVNROOT/Manuals/SMV_Verification_Guide/SMV_Verification_Guide.pdf
SMV_UG_GUIDE=$FDS_SVNROOT/Manuals/SMV_User_Guide/SMV_User_Guide.pdf
GEOM_NOTES=$FDS_SVNROOT/Manuals/FDS_User_Guide/geom_notes.pdf
GUIDE_DIR=$OUTPUT_DIR/Newest_Smokeview_Guides
UPDATE_GUIDES=$GUIDE_DIR/update_guides

THIS_FDS_AUTHOR=
THIS_FDS_FAILED=0
FDS_STATUS_FILE=$FDS_SVNROOT/FDS_status
LAST_FDS_FAILED=0
if [ -e $FDS_STATUS_FILE ] ; then
  LAST_FDS_FAILED=`cat $FDS_STATUS_FILE`
fi

# Load mailing list for status report
source $SMOKEBOT_DIR/firebot_email_list.sh

mailTo=$mailToSMV
if [[ "$LAST_FDS_FAILED" == "1" ]] ; then
  mailTo=$mailToFDS
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
  SVN_STATUSDIR=~/.smokebot
  SMV_SOURCE=$FDS_SVNROOT/SMV/source
  SVN_SMVFILE=$SVN_STATUSDIR/smv_revision
  SVN_SMVLOG=$SVN_STATUSDIR/smv_log

  FDS_SOURCE=$FDS_SVNROOT/FDS_Source
  SVN_FDSFILE=$SVN_STATUSDIR/fds_revision
  SVN_FDSLOG=$SVN_STATUSDIR/FDS_log

  MESSAGE_FILE=$SVN_STATUSDIR/message

  MKDIR $SVN_STATUSDIR
  cd $SMV_SOURCE
  svn update > /dev/null
  THIS_SMVSVN=`svn info | tail -3 | head -1 | awk '{print $4}'`
  THIS_SMVAUTHOR=`svn info | tail -4 | head -1 | awk '{print $4}'`
  LAST_SMVSVN=`cat $SVN_SMVFILE`
  svn log -r $THIS_SMVSVN > $SVN_SMVLOG

  cd $FDS_SOURCE
  svn update > /dev/null
  THIS_FDSSVN=`svn info | tail -3 | head -1 | awk '{print $4}'`
  THIS_FDSAUTHOR=`svn info | tail -4 | head -1 | awk '{print $4}'`
  LAST_FDSSVN=`cat $SVN_FDSFILE`
  svn log -r $THIS_FDSSVN > $SVN_FDSLOG

  if [[ $THIS_SMVSVN == $LAST_SMVSVN && $THIS_FDSSVN == $LAST_FDSSVN ]] ; then
    exit
  fi

  rm -f $MESSAGE_FILE
  if [[ $THIS_SMVSVN != $LAST_SMVSVN ]] ; then
    echo $THIS_SMVSVN>$SVN_SMVFILE
    echo -e "smokeview source has changed. $LAST_SMVSVN->$THIS_SMVSVN($THIS_SMVAUTHOR)" >> $MESSAGE_FILE
    cat $SVN_SMVLOG >> $MESSAGE_FILE
  fi
  if [[ $THIS_FDSSVN != $LAST_FDSSVN ]] ; then
    echo $THIS_FDSSVN>$SVN_FDSFILE
    echo -e "FDS source has changed. $LAST_FDSSVN->$THIS_FDSSVN($THIS_FDSAUTHOR)" >> $MESSAGE_FILE
    cat $SVN_FDSLOG >> $MESSAGE_FILE
#    RUNFDSCASES=1
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
   cd $FDS_SVNROOT
   chmod -R go+r *
}

clean_smokebot_history()
{
   
   # Clean Smokebot metafiles
   MKDIR $SMOKEBOT_DIR
   cd $SMOKEBOT_DIR
   MKDIR guides
   MKDIR history
   MKDIR output
   rm -rf output/* > /dev/null
   MKDIR $GUIDE_DIR
}

delete_unversioned_files()
{
   # Delete all unversioned SVN files
   svn status --no-ignore | grep '^[I?]' | cut -c 9- | while IFS= read -r f; do rm -rf "$f"; done
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
   if [ -e "$CFAST_SVNROOT" ]
   # If yes, then update the CFAST repository and compile CFAST
   then
      echo "Updating and compiling CFAST:" > $OUTPUT_DIR/stage0_cfast
      cd $CFAST_SVNROOT/CFAST
      svn revert -Rq *
      delete_unversioned_files

      # Update to latest SVN revision
      svn update >> $OUTPUT_DIR/stage0_cfast 2>&1
      
   # If no, then checkout the CFAST repository and compile CFAST
   else
      echo "Downloading and compiling CFAST:" > $OUTPUT_DIR/stage0_cfast
      mkdir -p $CFAST_SVNROOT
      cd $CFAST_SVNROOT

      svn co https://cfast.googlecode.com/svn/trunk/cfast/trunk/CFAST CFAST >> $OUTPUT_DIR/stage0_cfast 2>&1
      
   fi
    # Build CFAST
    cd $CFAST_SVNROOT/CFAST/intel_${platform}_64
    rm -f cfast6_${platform}_64
    make --makefile ../makefile clean &> /dev/null
    ./make_cfast.sh >> $OUTPUT_DIR/stage0_cfast 2>&1

   # Check for errors in CFAST compilation
   cd $CFAST_SVNROOT/CFAST/intel_${platform}_64
   if [ -e "cfast6_${platform}_64" ]
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
#  = Stage 1 - SVN operations =
#  ============================

clean_svn_repo()
{
   # Check to see if FDS repository exists
   if [ -e "$FDS_SVNROOT" ]
   then
      cd $FDS_SVNROOT
      svn revert -Rq *
      delete_unversioned_files
   # If not, create FDS repository and checkout
     dummy=true
   else
      echo "Downloading FDS repository:" >> $OUTPUT_DIR/stage1 2>&1
      cd $SMOKEBOT_HOME_DIR
      svn co https://fds-smv.googlecode.com/svn/trunk/FDS/trunk/ FDS-SMV >> $OUTPUT_DIR/stage1 2>&1
   fi
}

do_svn_checkout()
{
   cd $FDS_SVNROOT
   echo "Checking out latest revision." >> $OUTPUT_DIR/stage1 2>&1
   svn update >> $OUTPUT_DIR/stage1 2>&1
   SVN_REVISION=`tail -n 1 $OUTPUT_DIR/stage1 | sed "s/[^0-9]//g"`
}

check_svn_checkout()
{
   cd $FDS_SVNROOT
   # Check for SVN errors
   if [[ `grep -E 'Updated|At revision' $OUTPUT_DIR/stage1 | wc -l` -ne 1 ]];
   then
      echo "Errors from Stage 1 - SVN operations:" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage1 >> $ERROR_LOG
      echo "" >> $ERROR_LOG
      email_build_status
      exit
   else
      stage1_success=true
   fi
}

#  ==================================
#  = Stage 2a/b - Compile FDS debug =
#  ==================================

compile_fds_db()
{
   # Clean and compile FDS debug
   cd $FDS_SVNROOT/FDS_Compilation/${OPENMP}intel_${platform}_64_db
   rm -f fds_${OPENMP}intel_${platform}_64_db
   make --makefile ../makefile clean &> /dev/null
   ./make_fds.sh &> $OUTPUT_DIR/stage2a
}

compile_fds_mpi_db()
{
   # Clean and compile mpi FDS debug
   cd $FDS_SVNROOT/FDS_Compilation/mpi_intel_${platform}_64$IB$DB
   rm -f fds_mpi_intel_${platform}_64$IB$DB
   make --makefile ../makefile clean &> /dev/null
   ./make_fds.sh &> $OUTPUT_DIR/stage2b
}

check_compile_fds_db()
{
   # Check for errors in FDS debug compilation
   cd $FDS_SVNROOT/FDS_Compilation/${OPENMP}intel_${platform}_64_db
   if [ -e "fds_${OPENMP}intel_${platform}_64_db" ]
   then
      stage2a_success=true
   else
      echo "Errors from Stage 2a - Compile FDS debug:" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage2a >> $ERROR_LOG
      echo "" >> $ERROR_LOG
      THIS_FDS_FAILED=1
   fi

   # Check for compiler warnings/remarks
   if [[ `grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage2a` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 2a warnings:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage2a >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   # if the executable does not exist then an email has already been sent
      if [ -e "fds_${OPENMP}intel_${platform}_64_db" ] ; then
        THIS_FDS_FAILED=1
      fi
   fi
}

check_compile_fds_mpi_db()
{
   # Check for errors in FDS debug compilation
   cd $FDS_SVNROOT/FDS_Compilation/mpi_intel_${platform}_64$IB$DB
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
   cd $FDS_SVNROOT/Verification
   find .                        -name '*.stop' -exec rm -f {} \;
   find .                        -name '*.err' -exec rm -f {} \;
   find scripts/Outfiles         -name '*.out' -exec rm -f {} \;
   find Visualization            -name '*.smv' -exec rm -f {} \;
   find Immersed_Boundary_Method -name '*.smv' -exec rm -f {} \;
   find WUI                      -name '*.smv' -exec rm -f {} \;

   #  =====================
   #  = Run all SMV cases =
   #  =====================

   cd $FDS_SVNROOT/Verification/scripts

   # Submit SMV verification cases and wait for them to start
   echo 'Running SMV verification cases:' >> $OUTPUT_DIR/stage3a 2>&1
   ./Run_SMV_Cases.sh $USEINSTALL2 -m 2 -d -q $SMOKEBOT_QUEUE -j $JOBPREFIX >> $OUTPUT_DIR/stage3a 2>&1
#   ./Run_SMV_Cases.sh -S $USEINSTALL2 -m 2 -d -q $SMOKEBOT_QUEUE -j $JOBPREFIX >> $OUTPUT_DIR/stage3a 2>&1
#   ./Run_SMV_Cases.sh -M $USEINSTALL2 -m 2 -d -q $SMOKEBOT_QUEUE -j $JOBPREFIX >> $OUTPUT_DIR/stage3a 2>&1

   # Wait for SMV verification cases to end
   wait_verification_cases_debug_end

}

check_verification_cases_debug()
{
   # Scan and report any errors in FDS verification cases
   cd $FDS_SVNROOT/Verification

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

#  ======================================================
#  = Stage 3b - Run FDS verification cases (debug mode) =
#  ======================================================

wait_fds_verification_cases_debug_end()
{
   # Scans qstat and waits for verification cases to end
   if [[ "$SMOKEBOT_QUEUE" == "none" ]]
   then
     while [[ `ps -u $USER -f | fgrep .fds | grep -v grep` != '' ]]; do
        JOBS_REMAINING=`ps -u $USER -f | fgrep .fds | grep -v grep | wc -l`
        echo "Waiting for ${JOBS_REMAINING} verification cases to complete." >> $OUTPUT_DIR/stage3b
        TIME_LIMIT_STAGE="3"
        check_time_limit
        sleep 30
     done
   else
     while [[ `qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX` != '' ]]; do
        JOBS_REMAINING=`qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX | wc -l`
        echo "Waiting for ${JOBS_REMAINING} verification cases to complete." >> $OUTPUT_DIR/stage3b
        TIME_LIMIT_STAGE="3"
        check_time_limit
        sleep 30
     done
   fi
}

run_fds_verification_cases_debug()
{
   #  ======================
   #  = Remove .stop files =
   #  ======================

   # Remove all .stop and .err files from Verification directories (recursively)
   cd $FDS_SVNROOT/Verification
   find . -name '*.stop' -exec rm -f {} \;
   find . -name '*.err' -exec rm -f {} \;
   find . -name '*.out' -exec rm -f {} \;
   find . -name '*.smv' -exec rm -f {} \;
   find . -name '*.smv' -exec rm -f {} \;
   find . -name '*.smv' -exec rm -f {} \;

   #  =====================
   #  = Run all FDS cases =
   #  =====================

   cd $FDS_SVNROOT/Verification

   # Submit FDS verification cases and wait for them to start
   echo 'Running FDS verification cases:' >> $OUTPUT_DIR/stage3b 2>&1
   ./Run_FDS_Cases.sh -m 2 -d -q $SMOKEBOT_QUEUE -j $JOBPREFIX >> $OUTPUT_DIR/stage3b 2>&1

   # Wait for FDS verification cases to end
   wait_fds_verification_cases_debug_end

}

check_fds_verification_cases_debug()
{
   # Scan and report any errors in FDS verification cases
   cd $FDS_SVNROOT/Verification

   if [[ `grep -rIi 'Run aborted' $OUTPUT_DIR/stage3b` == "" ]] && \
      [[ `grep -rIi 'Segmentation' Visualization/* WUI/* Immersed_Boundary_Method/*` == "" ]] && \
      [[ `grep -rI 'ERROR:' Visualization/* WUI/* Immersed_Boundary_Method/*` == "" ]] && \
      [[ `grep -rIi 'STOP: Numerical' Visualization/* WUI/* Immersed_Boundary_Method/*` == "" ]] && \
      [[ `grep -rIi -A 20 'forrtl' Visualization/* WUI/* Immersed_Boundary_Method/*` == "" ]]
   then
      stage3b_success=true
   else
      grep -rIi 'Run aborted' $OUTPUT_DIR/stage3b > $OUTPUT_DIR/stage3b_errors
      grep -rIi 'Segmentation' Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage3b_errors
      grep -rI 'ERROR:' Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage3b_errors
      grep -rIi 'STOP: Numerical' -rIi Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage3b_errors
      grep -rIi -A 20 'forrtl' Visualization/* WUI/* Immersed_Boundary_Method/* >> $OUTPUT_DIR/stage3b_errors
      
      echo "Errors from Stage 3b - Run verification cases (debug mode):" >> $ERROR_LOG
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

#  ====================================
#  = Stage 4a/b - Compile FDS release =
#  ====================================

compile_fds()
{
   # Clean and compile FDS
   cd $FDS_SVNROOT/FDS_Compilation/${OPENMP}intel_${platform}_64
   rm -f fds_${OPENMP}intel_${platform}_64
   make --makefile ../makefile clean &> /dev/null
   ./make_fds.sh &> $OUTPUT_DIR/stage4a
}

compile_fds_mpi()
{
   # Clean and compile FDS
   cd $FDS_SVNROOT/FDS_Compilation/mpi_intel_${platform}_64$IB
   rm -f fds_mpi_intel_${platform}_64$IB
   make --makefile ../makefile clean &> /dev/null
   ./make_fds.sh &> $OUTPUT_DIR/stage4b
}

check_compile_fds()
{
   # Check for errors in FDS compilation
   cd $FDS_SVNROOT/FDS_Compilation/${OPENMP}intel_${platform}_64
   if [ -e "fds_${OPENMP}intel_${platform}_64" ]
   then
      stage4a_success=true
   else
      echo "Errors from Stage 4a - Compile FDS release:" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage4a >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   # 'performing multi-file optimizations' and 'generating object file' are part of a normal compile
   if [[ `grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage4a | grep -v 'performing multi-file optimizations' | grep -v 'generating object file'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 4a warnings:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' $OUTPUT_DIR/stage4a | grep -v 'performing multi-file optimizations' | grep -v 'generating object file'>> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

check_compile_fds_mpi()
{
   # Check for errors in FDS compilation
   cd $FDS_SVNROOT/FDS_Compilation/mpi_intel_${platform}_64$IB
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

   # smokeview libraries
   cd $FDS_SVNROOT/SMV/Build/LIBS/lib_${platform}_intel_64
   echo 'Building Smokeview libraries:' >> $OUTPUT_DIR/stage5pre 2>&1
   ./makelibs.sh >> $OUTPUT_DIR/stage5pre 2>&1

   # smokezip:
   cd $FDS_SVNROOT/Utilities/smokezip/intel_${platform}_64
   rm -f *.o smokezip_${platform}_64
   echo 'Compiling smokezip:' >> $OUTPUT_DIR/stage5pre 2>&1
   ./make_zip.sh >> $OUTPUT_DIR/stage5pre 2>&1
   echo "" >> $OUTPUT_DIR/stage5pre 2>&1
   
   # smokediff:
   cd $FDS_SVNROOT/Utilities/smokediff/intel_${platform}_64
   rm -f *.o smokediff_${platform}_64
   echo 'Compiling smokediff:' >> $OUTPUT_DIR/stage5pre 2>&1
   ./make_diff.sh >> $OUTPUT_DIR/stage5pre 2>&1
   echo "" >> $OUTPUT_DIR/stage5pre 2>&1
   
   # background:
   cd $FDS_SVNROOT/Utilities/background/intel_${platform}_32
   rm -f *.o background
   echo 'Compiling background:' >> $OUTPUT_DIR/stage5pre 2>&1
   ./make_background.sh >> $OUTPUT_DIR/stage5pre 2>&1
   
  # wind2fds:
   cd $FDS_SVNROOT/Utilities/wind2fds/intel_${platform}_64
   rm -f *.o wind2fds_${platform}_64
   echo 'Compiling wind2fds:' >> $OUTPUT_DIR/stage5pre 2>&1
   ./make_wind.sh >> $OUTPUT_DIR/stage5pre 2>&1
   echo "" >> $OUTPUT_DIR/stage5pre 2>&1
   else
   echo "Warning: smokeview and utilities not built - C compiler not available" >> $OUTPUT_DIR/stage5pre 2>&1
   fi
}

is_file_installed()
{
  program=$1
  notfound=`$program -help |& tail -1 |& grep "not found" | wc -l`
  if [ "$notfound" == "1" ] ; then
    stage5pre_success="0"
    echo "***error: $program not installed" >> $OUTPUT_DIR/stage5pre
  fi
}

check_smv_utilities()
{
   if [ "$haveCC" == "1" ] ; then
     # Check for errors in SMV utilities compilation
     cd $FDS_SVNROOT
     if [ -e "$FDS_SVNROOT/Utilities/smokezip/intel_${platform}_64/smokezip_${platform}_64" ]  && \
        [ -e "$FDS_SVNROOT/Utilities/smokediff/intel_${platform}_64/smokediff_${platform}_64" ]  && \
        [ -e "$FDS_SVNROOT/Utilities/wind2fds/intel_${platform}_64/wind2fds_${platform}_64" ]  && \
        [ -e "$FDS_SVNROOT/Utilities/background/intel_${platform}_32/background" ]
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
   cd $FDS_SVNROOT/Verification
   find .                        -name '*.stop' -exec rm -f {} \;
   find .                        -name '*.err' -exec rm -f {} \;
   find scripts/Outfiles         -name '*.out' -exec rm -f {} \;
   find Visualization            -name '*.smv' -exec rm -f {} \;
   find Immersed_Boundary_Method -name '*.smv' -exec rm -f {} \;
   find WUI                      -name '*.smv' -exec rm -f {} \;

   # Start running all SMV verification cases
   cd $FDS_SVNROOT/Verification/scripts
   echo 'Running SMV verification cases:' >> $OUTPUT_DIR/stage5 2>&1
   ./Run_SMV_Cases.sh $USEINSTALL2 $RUN_OPENMP -q $SMOKEBOT_QUEUE -j $JOBPREFIX >> $OUTPUT_DIR/stage5 2>&1
#   ./Run_SMV_Cases.sh -S $USEINSTALL2 $RUN_OPENMP -q $SMOKEBOT_QUEUE -j $JOBPREFIX >> $OUTPUT_DIR/stage5 2>&1
#   ./Run_SMV_Cases.sh -M $USEINSTALL2 $RUN_OPENMP -q $SMOKEBOT_QUEUE -j $JOBPREFIX >> $OUTPUT_DIR/stage5 2>&1

   # Wait for all verification cases to end
   wait_verification_cases_release_end
}

check_verification_cases_release()
{
   # Scan and report any errors in FDS verification cases
   cd $FDS_SVNROOT/Verification

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
   # Clean and compile SMV debug
   cd $FDS_SVNROOT/SMV/Build/intel_${platform}_64
   rm -f smokeview_${platform}_64_db
   ./make_smv_db.sh &> $OUTPUT_DIR/stage6a
   fi
}

check_compile_smv_db()
{
   if [ "$haveCC" == "1" ] ; then
   # Check for errors in SMV debug compilation
   cd $FDS_SVNROOT/SMV/Build/intel_${platform}_64
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
   cd $FDS_SVNROOT/Verification/scripts
   ./Make_SMV_Pictures.sh $USEINSTALL -d 2>&1 | grep -v FreeFontPath &> $OUTPUT_DIR/stage6b
}

check_smv_pictures_db()
{
   # Scan and report any errors in make SMV pictures process
   cd $SMOKEBOT_DIR
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
   cd $SMOKEBOT_DIR
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
   # Clean and compile SMV
   cd $FDS_SVNROOT/SMV/Build/intel_${platform}_64
   rm -f smokeview_${platform}_64
   ./make_smv.sh $TESTFLAG &> $OUTPUT_DIR/stage6c
   fi
}

check_compile_smv()
{
   if [ "$haveCC" == "1" ] ; then
   # Check for errors in SMV release compilation
   cd $FDS_SVNROOT/SMV/Build/intel_${platform}_64
   if [ -e "smokeview_${platform}_64" ]
   then
      stage6c_success=true
   else
      echo "Errors from Stage 6c - Compile SMV release:" >> $ERROR_LOG
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
   cd $FDS_SVNROOT/Verification/scripts
   ./Make_SMV_Pictures.sh $TESTFLAG $USEINSTALL 2>&1 | grep -v FreeFontPath &> $OUTPUT_DIR/stage6d
}

check_smv_pictures()
{
   # Scan and report any errors in make SMV pictures process
   cd $SMOKEBOT_DIR
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
   cd $FDS_SVNROOT/Verification
   scripts/Make_SMV_Movies.sh 2>&1  &> $OUTPUT_DIR/stage6e
}

check_smv_movies()
{
   cd $SMOKEBOT_DIR
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
   cd $SMOKEBOT_DIR
   if [[ `grep -I -E "Warning" $OUTPUT_DIR/stage6d` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 6d - Make SMV pictures (release mode):" >> $WARNING_LOG
      grep -I -E "Warning" $OUTPUT_DIR/stage6d >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi

}

#  ======================================
#  = Stage 7 - FDS run time statistics =
#  ======================================

generate_timing_stats()
{
   cd $FDS_SVNROOT/Verification/scripts/
   export QFDS="$FDS_SVNROOT/Verification/scripts/copyout.sh"
   export RUNCFAST="$FDS_SVNROOT/Verification/scripts/copyout.sh"
   export RUNTFDS="$FDS_SVNROOT/Verification/scripts/copyout.sh"

   cd $FDS_SVNROOT/Verification
   scripts/SMV_Cases.sh
   scripts/SMV_geom_Cases.sh

   cd $FDS_SVNROOT/Utilities/Scripts
   ./fds_timing_stats.sh smokebot
}

archive_timing_stats()
{
   cd $FDS_SVNROOT/Utilities/Scripts
   cp fds_timing_stats.csv "$SMOKEBOT_DIR/history/${SVN_REVISION}_timing.csv"
}

#  ==================================
#  = Stage 8 - Build FDS-SMV Guides =
#  ==================================

check_guide()
{
   # Scan and report any errors in build process for guides
   SMOKEBOT_MANDIR=/var/www/html/smokebot/manuals/
   cd $SMOKEBOT_DIR
   if [[ `grep "! LaTeX Error:" -I $1` == "" ]]
   then
      if [ -d $SMOKEBOT_MANDIR ] ; then
        cp $2 $SMOKEBOT_MANDIR/.
      fi
      if [ -d $SMV_Summary/manuals ] ; then
        cp $2 $SMV_Summary/manuals/.
      fi
      cp $2 $GUIDE_DIR/.
      chmod 664 $GUIDE_DIR/$2
   else
      echo "Errors from Stage 8 - Build FDS-SMV Guides:" >> $ERROR_LOG
      echo $3 >> $ERROR_LOG
      grep "! LaTeX Error:" -I $1 >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for LaTeX warnings (undefined references or duplicate labels)
   if [[ `grep -E "undefined|multiply defined|multiply-defined" -I ${1}` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 8 warnings:" >> $WARNING_LOG
      echo $3 >> $WARNING_LOG
      grep -E "undefined|multiply defined|multiply-defined" -I $1 >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

make_guide()
{
   document=$1
   directory=$2
   label=$3

   cd $directory
   export TEXINPUTS=".:../LaTeX_Style_Files:"
   pdflatex -interaction nonstopmode $document &> $OUTPUT_DIR/stage8_$document
   bibtex $document &> $OUTPUT_DIR/stage8_$document
   pdflatex -interaction nonstopmode $document &> $OUTPUT_DIR/stage8_$document
   pdflatex -interaction nonstopmode $document &> $OUTPUT_DIR/stage8_$document

   # Check guide for completion and copy to website if successful
   check_guide $OUTPUT_DIR/stage8_$document $directory/$document.pdf $label
}

#  =====================================================
#  = Build status reporting - email and save functions =
#  =====================================================

save_build_status()
{
   cd $SMOKEBOT_DIR
   # Save status outcome of build to a text file
   if [[ -e $WARNING_LOG && -e $ERROR_LOG ]]
   then
     cat "" >> $ERROR_LOG
     cat $WARNING_LOG >> $ERROR_LOG
     echo "Build failure and warnings for Revision ${SVN_REVISION}." > "$SMOKEBOT_DIR/history/${SVN_REVISION}.txt"
     cat $ERROR_LOG > "$SMOKEBOT_DIR/history/${SVN_REVISION}_errors.txt"
     touch output/status_errors_and_warnings

   # Check for errors only
   elif [ -e $ERROR_LOG ]
   then
      echo "Build failure for Revision ${SVN_REVISION}." > "$SMOKEBOT_DIR/history/${SVN_REVISION}.txt"
      cat $ERROR_LOG > "$SMOKEBOT_DIR/history/${SVN_REVISION}_errors.txt"
      touch output/status_errors

   # Check for warnings only
   elif [ -e $WARNING_LOG ]
   then
      echo "Revision ${SVN_REVISION} has warnings." > "$SMOKEBOT_DIR/history/${SVN_REVISION}.txt"
      cat $WARNING_LOG > "$SMOKEBOT_DIR/history/${SVN_REVISION}_warnings.txt"
      touch output/status_warnings

   # No errors or warnings
   else
      echo "Build success! Revision ${SVN_REVISION} passed all build tests." > "$SMOKEBOT_DIR/history/${SVN_REVISION}.txt"
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
   echo ".        total: $DIFF_SCRIPT_TIME" >> $TIME_LOG
   echo ".      results: http://$WEBHOSTNAME/VV/SMV2" >> $TIME_LOG
  if [[ $THIS_SMVSVN != $LAST_SMVSVN ]] ; then
    cat $SVN_SMVLOG >> $TIME_LOG
  fi
  if [[ $THIS_FDSSVN != $LAST_FDSSVN ]] ; then
    cat $SVN_FDSLOG >> $TIME_LOG
  fi
   echo "----------------------------------------------" >> $TIME_LOG
   cd $SMOKEBOT_DIR
   # Check for warnings and errors
   if [[ -e $WARNING_LOG && -e $ERROR_LOG ]]
   then
     cat $TIME_LOG >> $ERROR_LOG
     cat $TIME_LOG >> $WARNING_LOG
     # Send email with failure message and warnings, body of email contains appropriate log file
     mail -s "smokebot build failure and warnings on ${hostname}. Revision ${SVN_REVISION}." $mailTo < $ERROR_LOG > /dev/null

   # Check for errors only
   elif [ -e $ERROR_LOG ]
   then
     cat $TIME_LOG >> $ERROR_LOG
      # Send email with failure message, body of email contains error log file
      mail -s "smokebot build failure on ${hostname}. Revision ${SVN_REVISION}." $mailTo < $ERROR_LOG > /dev/null

   # Check for warnings only
   elif [ -e $WARNING_LOG ]
   then
     cat $TIME_LOG >> $WARNING_LOG
      # Send email with success message, include warnings
      mail -s "smokebot build success with warnings on ${hostname}. Revision ${SVN_REVISION}." $mailTo < $WARNING_LOG > /dev/null

   # No errors or warnings
   else
      # Send empty email with success message
      if [ "$MUTT_MISSING" == "1" ] ; then
        mail -s "smokebot build success on ${hostname}! Revision ${SVN_REVISION}." $mailTo < $TIME_LOG > /dev/null
      else
        mail -s "smokebot build success on ${hostname}! Revision ${SVN_REVISION}." $mailTo < $TIME_LOG > /dev/null
#        mutt -a $SMV_VG_GUIDE -a $SMV_UG_GUIDE -s "smokebot build success on ${hostname}! Revision ${SVN_REVISION}." $mailTo < $TIME_LOG > /dev/null
      fi
      touch $UPDATE_GUIDES
      chmod 644 $UPDATE_GUIDES
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
clean_svn_repo
do_svn_checkout
check_svn_checkout
PRELIM_end=`GET_TIME`
DIFF_PRELIM=`GET_DURATION $PRELIM_beg $PRELIM_end`
echo "Preliminary: $DIFF_PRELIM" >> $STAGE_STATUS

### Stage 2a ###
BUILDFDS_beg=`GET_TIME`
#compile_fds_db
#check_compile_fds_db
stage2a_success=true

### Stage 2b ###
compile_fds_mpi_db
check_compile_fds_mpi_db

### Stage 4a ###
stage4_beg=`GET_TIME`
#if [[ $stage2a_success ]] ; then
#   compile_fds
#   check_compile_fds
#fi
stage4a_success=true

### Stage 4b ###
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
if [[ $stage2a_success && $stage2b_success && "$RUNDEBUG" == "1" ]] ; then
   run_verification_cases_debug
   check_verification_cases_debug
   if [ "$RUNFDSCASES" == "1" ] ; then
     run_fds_verification_cases_debug
     check_fds_verification_cases_debug
   fi
fi

### Stage 5 ###
if [[ $stage4a_success && $stage4b_success ]] ; then
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

### Stage 6b ###
MAKEPICTURES_beg=`GET_TIME`
if [[ $stage4a_success && $stage4b_success && $stage6a_success && "$RUNDEBUG" == "1" ]] ; then
# for now skip image generation in debug mode
DUMMY=
#  make_smv_pictures_db
#  check_smv_pictures_db
fi

### Stage 6d ###
if [[ $stage4a_success && $stage4b_success && $stage6c_success ]] ; then
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
if [[ $stage4a_success && $stage4b_success ]] ; then
  generate_timing_stats
  archive_timing_stats
fi

### Stage 8 ###
MAKEGUIDES_beg=`GET_TIME`
if [[ $stage4a_success && $stage4b_success && $stage6d_success ]] ; then
#  make_guide geom_notes $FDS_SVNROOT/Manuals/FDS_User_Guide 'geometry notes'
  make_guide SMV_User_Guide $FDS_SVNROOT/Manuals/SMV_User_Guide 'SMV User Guide'
  make_guide SMV_Technical_Reference_Guide $FDS_SVNROOT/Manuals/SMV_Technical_Reference_Guide 'SMV Technical Reference Guide'
  make_guide SMV_Verification_Guide $FDS_SVNROOT/Manuals/SMV_Verification_Guide 'SMV Verification Guide'
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
