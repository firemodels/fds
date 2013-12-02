#!/bin/bash

# Firebot
# FDS Automatic Verification and Validation Test Bot
# Kristopher Overholt
# 6/22/2012

# The Firebot script is part of an automated continuous integration system.
# Please consult the Utilities/Firebot/README.txt file and the
# FDS Configuration Management Plan for more information.

#  ===================
#  = Input variables =
#  ===================

# Firebot mode (verification or validation); default mode: verification
FIREBOT_MODE="verification"

# Firebot's username
FIREBOT_USERNAME="firebot"

# Change to home directory
cd
FIREBOT_HOME_DIR="`pwd`"

# Set unlimited stack size
ulimit -s unlimited

# Additional definitions
FIREBOT_DIR="$FIREBOT_HOME_DIR/firebot"
FDS_SVNROOT="$FIREBOT_HOME_DIR/FDS-SMV"
CFAST_SVNROOT="$FIREBOT_HOME_DIR/cfast"
TIME_LOG=$FIREBOT_DIR/output/timings
ERROR_LOG=$FIREBOT_DIR/output/errors
WARNING_LOG=$FIREBOT_DIR/output/warnings
DB=_db
IB=
if [ "$FDSNETWORK" == "infiniband" ] ; then
IB=ib
fi

# Load mailing list for status report
source $FIREBOT_DIR/firebot_email_list.sh

function usage {
echo "firebot.sh [ -q queue_name -r revision_number -s -u svn_username -y ]"
echo "Runs Firebot V&V testing script"
echo ""
echo "Options"
echo "-q - queue_name - run cases using the queue queue_name"
echo "     default: firebot"
echo ""
echo "-r - revision_number - run cases using a specific SVN revision number"
echo "     default: (none, latest SVN HEAD)"
echo ""
echo "-s - skip fixing SVN properties"
echo "     default: SKIP_SVN_PROPS is undefined (false)"
echo ""
echo "-u - specify SVN username to use"
echo "     default: fds.firebot"
echo ""
echo "-v - run Firebot in validation mode"
echo "     default: (none, Firebot runs in verification mode)"
echo ""
echo "-y - run Firebot as any user (warning!)"
echo "     default: (none)"
exit
}

QUEUE=firebot
SVN_REVISION=
SVN_USERNAME=fds.firebot
while getopts 'hq:r:su:vy' OPTION
do
case $OPTION in
  h)
   usage;
   ;;
  q)
   QUEUE="$OPTARG"
   ;;
  r)
   SVN_REVISION="$OPTARG"
   ;;
  s)
   SKIP_SVN_PROPS=true
   ;;
  u)
   SVN_USERNAME="$OPTARG"
   ;;
  v)
   FIREBOT_MODE="validation"
   # Disable time limit email
   TIME_LIMIT_EMAIL_NOTIFICATION="sent"
   mailToFDS="koverholt@gmail.com" # For validationbot testing
   ;;
  y)
   RUN_AS_ANOTHER_USER=true
   ;;
esac
done
shift $(($OPTIND-1))

#  ====================
#  = End user warning =
#  ====================

if [[ $RUN_AS_ANOTHER_USER ]] ; then
# Continue along
:
else
   # Warn if running as user other than firebot
   if [[ `whoami` == "$FIREBOT_USERNAME" ]];
      then
         # Continue along
         :
      else
         echo "Warning: You are running the Firebot script as an end user."
         echo "This script will definitely modify and/or erase your repository."
         echo "If you wish to continue, run Firebot with the -y option at your own risk."
         echo "Terminating script."
         exit
   fi
fi

#  =============================================
#  = Firebot timing and notification mechanism =
#  =============================================

# This routine checks the elapsed time of Firebot.
# If Firebot runs more than 12 hours, an email notification is sent.
# This is a notification only and does not terminate Firebot.
# This check runs during Stages 3 and 5.

# Start firebot timer
START_TIME=$(date +%s)

# Set time limit (43,200 seconds = 12 hours)
TIME_LIMIT=43200
TIME_LIMIT_EMAIL_NOTIFICATION="unsent"

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
         echo -e "Firebot has been running for more than 12 hours in Stage ${TIME_LIMIT_STAGE}. \n\nPlease ensure that there are no problems. \n\nThis is a notification only and does not terminate Firebot." | mail -s "[Firebot@$hostname] Notice: Firebot has been running for more than 12 hours." $mailToFDS > /dev/null
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

clean_firebot_history()
{
   # Clean Firebot metafiles
   cd $FIREBOT_DIR
   rm output/* > /dev/null
}

#  ========================
#  ========================
#  = Firebot Build Stages =
#  ========================
#  ========================

#  ===================================
#  = Stage 0 - External dependencies =
#  ===================================

update_and_compile_cfast()
{
   cd $FIREBOT_HOME_DIR

   # Check to see if CFAST repository exists
   if [ -e "$CFAST_SVNROOT" ]
   # If yes, then update the CFAST repository and compile CFAST
   then
      echo "Updating and compiling CFAST:" >> $FIREBOT_DIR/output/stage0_cfast
      cd $CFAST_SVNROOT/CFAST
      
      # Clean unversioned and modified files
      svn revert -Rq *
      svn status --no-ignore | grep '^[I?]' | cut -c 9- | while IFS= read -r f; do rm -rf "$f"; done
      
      # Update to latest SVN revision
      svn update >> $FIREBOT_DIR/output/stage0_cfast 2>&1
      
      # Build CFAST
      cd $CFAST_SVNROOT/CFAST/intel_linux_64
      make -f ../makefile clean &> /dev/null
      ./make_cfast.sh >> $FIREBOT_DIR/output/stage0_cfast 2>&1
   # If no, then checkout the CFAST repository and compile CFAST
   else
      echo "Downloading and compiling CFAST:" >> $FIREBOT_DIR/output/stage0_cfast
      mkdir -p $CFAST_SVNROOT
      cd $CFAST_SVNROOT

      # Checkout latest CFAST SVN revision
      svn co https://cfast.googlecode.com/svn/trunk/cfast/trunk/CFAST CFAST >> $FIREBOT_DIR/output/stage0_cfast 2>&1
      
      # Build CFAST
      cd $CFAST_SVNROOT/CFAST/intel_linux_64
      make -f ../makefile clean &> /dev/null
      ./make_cfast.sh >> $FIREBOT_DIR/output/stage0_cfast 2>&1
   fi

   # Check for errors in CFAST compilation
   cd $CFAST_SVNROOT/CFAST/intel_linux_64
   if [ -e "cfast6_linux_64" ]
   then
      stage0_success=true
   else
      echo "Errors from Stage 0 - CFAST:" >> $ERROR_LOG
      echo "CFAST failed to compile" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage0_cfast >> $ERROR_LOG
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
   # If yes, clean FDS repository
   then
      # Revert and clean up temporary unversioned and modified versioned repository files
      cd $FDS_SVNROOT
      svn revert -Rq *
      svn status --no-ignore | grep '^[I?]' | cut -c 9- | while IFS= read -r f; do rm -rf "$f"; done
   # If not, create FDS repository and checkout
   else
      echo "Downloading FDS repository:" >> $FIREBOT_DIR/output/stage1 2>&1
      cd $FIREBOT_HOME_DIR
      svn co https://fds-smv.googlecode.com/svn/trunk/FDS/trunk/ FDS-SMV --username $SVN_USERNAME >> $FIREBOT_DIR/output/stage1 2>&1
   fi
}

do_svn_checkout()
{
   cd $FDS_SVNROOT
   # If an SVN revision number is specified, then get that revision
   if [[ $SVN_REVISION != "" ]]; then
      echo "Checking out revision r${SVN_REVISION}." >> $FIREBOT_DIR/output/stage1 2>&1
      svn update -r $SVN_REVISION >> $FIREBOT_DIR/output/stage1 2>&1
   # If no SVN revision number is specified, then get the latest revision
   else
      echo "Checking out latest revision." >> $FIREBOT_DIR/output/stage1 2>&1
      svn update >> $FIREBOT_DIR/output/stage1 2>&1

      # Only run if firebot is in "verification" mode
      if [ $FIREBOT_MODE == "verification" ] ; then
         # Bump SVN revision number for all guides (so that the SVN revision keyword gets updated)
         echo "Bump SVN revision number for all guides." >> $FIREBOT_DIR/output/stage1 2>&1
         CURRENT_TIMESTAMP=`date`
         sed -i "s/.*% dummy comment to force svn change.*/% dummy comment to force svn change - ${CURRENT_TIMESTAMP}/" $FDS_SVNROOT/Manuals/FDS_User_Guide/FDS_User_Guide.tex
         sed -i "s/.*% dummy comment to force svn change.*/% dummy comment to force svn change - ${CURRENT_TIMESTAMP}/" $FDS_SVNROOT/Manuals/FDS_Technical_Reference_Guide/FDS_Technical_Reference_Guide.tex
         sed -i "s/.*% dummy comment to force svn change.*/% dummy comment to force svn change - ${CURRENT_TIMESTAMP}/" $FDS_SVNROOT/Manuals/FDS_Verification_Guide/FDS_Verification_Guide.tex
         sed -i "s/.*% dummy comment to force svn change.*/% dummy comment to force svn change - ${CURRENT_TIMESTAMP}/" $FDS_SVNROOT/Manuals/FDS_Validation_Guide/FDS_Validation_Guide.tex
         sed -i "s/.*% dummy comment to force svn change.*/% dummy comment to force svn change - ${CURRENT_TIMESTAMP}/" $FDS_SVNROOT/Manuals/FDS_Configuration_Management_Plan/FDS_Configuration_Management_Plan.tex

         # Commit back results
         svn commit -m 'Firebot: Bump SVN revision number for all guides' &> /dev/null
      fi
      
      echo "Re-checking out latest revision." >> $FIREBOT_DIR/output/stage1 2>&1
      svn update >> $FIREBOT_DIR/output/stage1 2>&1
      SVN_REVISION=`tail -n 1 $FIREBOT_DIR/output/stage1 | sed "s/[^0-9]//g"`
   fi
}

check_svn_checkout()
{
   cd $FDS_SVNROOT
   # Check for SVN errors
   if [[ `grep -E 'Updated|At revision' $FIREBOT_DIR/output/stage1 | wc -l` -ne 2 ]];
   then
      echo "Errors from Stage 1 - SVN operations:" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage1 >> $ERROR_LOG
      echo "" >> $ERROR_LOG
      if [ $FIREBOT_MODE == "verification" ] ; then
         email_build_status 'Firebot' 'Build'
      elif [ $FIREBOT_MODE == "validation" ] ; then
         email_build_status 'Validationbot' 'Validation'
      fi
      exit
   else
      stage1_success=true
   fi
}

fix_svn_properties()
{
   # This function fixes SVN properties
   # (e.g., svn:executable, svn:keywords, svn:eol-style, and svn:mime-type)
   # throughout the FDS-SMV repository.

   # cd to SVN root
   cd $FDS_SVNROOT

   # Delete all svn:executable properties
   svn propdel svn:executable --recursive &> /dev/null

   # Restore local executable property to svn-fix-props.pl
   chmod +x Utilities/Subversion/svn-fix-props.pl &> /dev/null

   # Run svn-fix-props.pl script (fixes all SVN properties)
   Utilities/Subversion/svn-fix-props.pl ./ &> /dev/null

   # Commit back results
   svn commit -m 'Firebot: Fix SVN properties throughout repository' &> /dev/null
}

#  ================================
#  = Stage 2a - Compile FDS debug =
#  ================================

compile_fds_db()
{
   # Clean and compile FDS debug
   cd $FDS_SVNROOT/FDS_Compilation/intel_linux_64_db
   make -f ../makefile clean &> /dev/null
   ./make_fds.sh &> $FIREBOT_DIR/output/stage2a
}

check_compile_fds_db()
{
   # Check for errors in FDS debug compilation
   cd $FDS_SVNROOT/FDS_Compilation/intel_linux_64_db
   if [ -e "fds_intel_linux_64_db" ]
   then
      stage2a_success=true
   else
      echo "Errors from Stage 2a - Compile FDS debug:" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage2a >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2a` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 2a - Compile FDS debug:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2a >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

#  ====================================
#  = Stage 2b - Compile FDS MPI debug =
#  ====================================

compile_fds_mpi_db()
{
   # Clean and compile FDS MPI debug
   cd $FDS_SVNROOT/FDS_Compilation/mpi_intel_linux_64$IB$DB
   make -f ../makefile clean &> /dev/null
   ./make_fds.sh &> $FIREBOT_DIR/output/stage2b
}

check_compile_fds_mpi_db()
{
   # Check for errors in FDS MPI debug compilation
   cd $FDS_SVNROOT/FDS_Compilation/mpi_intel_linux_64$IB$DB
   if [ -e "fds_mpi_intel_linux_64$IB$DB" ]
   then
      stage2b_success=true
   else
      echo "Errors from Stage 2b - Compile FDS MPI debug:" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage2b >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   # grep -v 'feupdateenv ...' ignores a known FDS MPI compiler warning (http://software.intel.com/en-us/forums/showthread.php?t=62806)
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2b | grep -v 'feupdateenv is not implemented'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 2b - Compile FDS MPI debug:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2b | grep -v 'feupdateenv is not implemented' >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

#  ===================================================
#  = Stage 2c - Compile and inspect FDS OpenMP debug =
#  ===================================================

compile_fds_openmp_db()
{
   # Clean and compile FDS OpenMP debug
   cd $FDS_SVNROOT/FDS_Compilation/openmp_intel_linux_64_db
   make -f ../makefile clean &> /dev/null
   ./make_fds.sh &> $FIREBOT_DIR/output/stage2c
}

check_compile_fds_openmp_db()
{
   # Check for errors in FDS OpenMP debug compilation
   cd $FDS_SVNROOT/FDS_Compilation/openmp_intel_linux_64_db
   if [ -e "fds_openmp_intel_linux_64_db" ]
   then
      stage2c_success=true
   else
      echo "Errors from Stage 2c - Compile and inspect FDS OpenMP debug:" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage2c >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2c` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 2c - Compile and inspect FDS OpenMP debug:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2c >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

inspect_fds_openmp_db()
{
   # Perform OpenMP thread checking (locate deadlocks and data races)
   cd $FDS_SVNROOT/Utilities/Scripts
   ./inspect_openmp.sh &> $FIREBOT_DIR/output/stage2c_inspect
}

check_inspect_fds_openmp_db()
{
   # Scan for errors in thread checking results
   cd $FDS_SVNROOT/Utilities/Scripts
   if [[ `grep -i -E 'problem' ${FIREBOT_DIR}/output/stage2c_inspect` == "" ]]
   then
      # Continue along
      :
   else
      echo "Errors from Stage 2c - Compile and inspect FDS OpenMP debug:" >> $ERROR_LOG
      cat ${FIREBOT_DIR}/output/stage2c_inspect >> $ERROR_LOG
      echo "" >> $ERROR_LOG
      echo "For more details, view the inspector log in the FDS-SMV/Utilities/Scripts folder" >> $ERROR_LOG
      echo "by using the FDS-SMV/Utilities/Scripts/inspect_report.sh script." >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi
}

#  =================================================
#  = Stage 3 - Run verification cases (debug mode) =
#  =================================================

select_validation_set()
{
   cd $FDS_SVNROOT/Validation

   # List five oldest validation sets in the $FDS_SVNROOT/Validation/Process_All_Output.sh script
   # based on the modification date of $VDIR/FDS_Output_Files
   # The result is an array of five validation sets from 0 to 4
   VALIDATION_SETS=(`grep '$VDIR' Process_All_Output.sh | grep -v "#" | xargs -n 1 dirname | xargs -n 1 dirname | xargs -n 1 basename | xargs -i ls --time-style full-iso -tl -d -1 {}"/FDS_Output_Files" | sort -k 6 | cut -d ' ' -f9 | head -n 5 | xargs -n 1 dirname`)

   # Select a random validation set number from 0 to 4
   RANDOM_VALIDATION_SET=$(($RANDOM % 5))
   CURRENT_VALIDATION_SET=${VALIDATION_SETS[RANDOM_VALIDATION_SET]}

   echo "Current Validation Set: ${CURRENT_VALIDATION_SET}" >> $FIREBOT_DIR/output/stage3
   echo "" >> $FIREBOT_DIR/output/stage3 2>&1
}

wait_cases_debug_start()
{
   # Scans qstat and waits for verification cases to start
   while [[ `qstat | grep $(whoami) | awk '{print $5}' | grep Q` != '' ]]; do
      JOBS_REMAINING=`qstat | grep $(whoami) | awk '{print $5}' | grep Q | wc -l`
      echo "Waiting for ${JOBS_REMAINING} ${1} cases to start." >> $FIREBOT_DIR/output/stage3
      TIME_LIMIT_STAGE="3"
      check_time_limit
      sleep 30
   done
}

wait_cases_debug_end()
{
   # Scans qstat and waits for verification cases to end
   while [[ `qstat | grep $(whoami)` != '' ]]; do
      JOBS_REMAINING=`qstat | grep $(whoami) | wc -l`
      echo "Waiting for ${JOBS_REMAINING} ${1} cases to complete." >> $FIREBOT_DIR/output/stage3
      TIME_LIMIT_STAGE="3"
      check_time_limit
      sleep 30
   done
}

run_verification_cases_debug()
{
   #  ============================
   #  = Run all FDS serial cases =
   #  ============================

   cd $FDS_SVNROOT/Verification

   # Submit FDS verification cases and wait for them to start (run serial cases in debug mode on firebot queue)
   echo 'Running FDS verification cases (serial):' >> $FIREBOT_DIR/output/stage3
   ./Run_FDS_Cases.sh -c serial -d -q $QUEUE >> $FIREBOT_DIR/output/stage3 2>&1
   wait_cases_debug_start 'verification'

   # Wait some additional time for all cases to start
   sleep 30

   # Stop all cases
   ./Run_FDS_Cases.sh -c serial -d -s >> $FIREBOT_DIR/output/stage3 2>&1
   echo "" >> $FIREBOT_DIR/output/stage3 2>&1

   # Wait for serial verification cases to end
   wait_cases_debug_end 'verification'

   #  =========================
   #  = Run all FDS MPI cases =
   #  =========================

   cd $FDS_SVNROOT/Verification

   # Submit FDS verification cases and wait for them to start (run MPI cases in debug mode on firebot queue)
   echo 'Running FDS verification cases (MPI):' >> $FIREBOT_DIR/output/stage3 2>&1
   ./Run_FDS_Cases.sh -c mpi -d -q $QUEUE >> $FIREBOT_DIR/output/stage3 2>&1
   wait_cases_debug_start 'verification'

   # Wait some additional time for all cases to start
   sleep 30

   # Stop all cases
   ./Run_FDS_Cases.sh -c mpi -d -s >> $FIREBOT_DIR/output/stage3 2>&1
   echo "" >> $FIREBOT_DIR/output/stage3 2>&1

   # Wait for MPI verification cases to end
   wait_cases_debug_end 'verification'

   #  =====================
   #  = Run all SMV cases =
   #  =====================

   cd $FDS_SVNROOT/Verification/scripts

   # Submit SMV verification cases and wait for them to start (run SMV cases in debug mode on firebot queue)
   echo 'Running SMV verification cases:' >> $FIREBOT_DIR/output/stage3 2>&1
   ./Run_SMV_Cases.sh -d -q $QUEUE >> $FIREBOT_DIR/output/stage3 2>&1
   wait_cases_debug_start 'verification'

   # Wait some additional time for all cases to start
   sleep 30

   # Stop all cases
   ./Run_SMV_Cases.sh -d -s >> $FIREBOT_DIR/output/stage3 2>&1
   echo "" >> $FIREBOT_DIR/output/stage3 2>&1

   # Wait for SMV verification cases to end
   wait_cases_debug_end 'verification'

   #  ======================
   #  = Remove .stop files =
   #  ======================

   # Remove all .stop files from Verification directories (recursively)
   cd $FDS_SVNROOT/Verification
   find . -name '*.stop' -exec rm -f {} \;
}

run_validation_cases_debug()
{
   # Placeholder
   sleep 30
}

check_cases_debug()
{
   # Scan and report any errors in FDS cases
   cd $1

   if [[ `grep -rI 'Run aborted' ${FIREBOT_DIR}/output/stage3` == "" ]] && \
      [[ `grep -rI Segmentation *` == "" ]] && \
      [[ `grep -rI ERROR: *` == "" ]] && \
      [[ `grep -rI 'STOP: Numerical' *` == "" ]] && \
      [[ `grep -rI -A 20 forrtl *` == "" ]]
   then
      stage3_success=true
   else
      grep -rI 'Run aborted' $FIREBOT_DIR/output/stage3 >> $FIREBOT_DIR/output/stage3_errors
      grep -rI Segmentation * >> $FIREBOT_DIR/output/stage3_errors
      grep -rI ERROR: * >> $FIREBOT_DIR/output/stage3_errors
      grep -rI 'STOP: Numerical' * >> $FIREBOT_DIR/output/stage3_errors
      grep -rI -A 20 forrtl * >> $FIREBOT_DIR/output/stage3_errors
      
      echo "Errors from Stage 3 - Run ${2} cases (debug mode):" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage3_errors >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # After Stage 3, delete all unversioned FDS output files before continuing
   cd $1
   svn status --no-ignore | grep '^[I?]' | cut -c 9- | while IFS= read -r f; do rm -rf "$f"; done
}

#  ==================================
#  = Stage 4a - Compile FDS release =
#  ==================================

compile_fds()
{
   # Clean and compile FDS
   cd $FDS_SVNROOT/FDS_Compilation/intel_linux_64
   make -f ../makefile clean &> /dev/null
   ./make_fds.sh &> $FIREBOT_DIR/output/stage4a
}

check_compile_fds()
{
   # Check for errors in FDS compilation
   cd $FDS_SVNROOT/FDS_Compilation/intel_linux_64
   if [ -e "fds_intel_linux_64" ]
   then
      stage4a_success=true
   else
      echo "Errors from Stage 4a - Compile FDS release:" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage4a >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   # 'performing multi-file optimizations' and 'generating object file' are part of a normal compile
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage4a | grep -v 'performing multi-file optimizations' | grep -v 'generating object file'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 4a - Compile FDS release:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage4a | grep -v 'performing multi-file optimizations' | grep -v 'generating object file' >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

#  ======================================
#  = Stage 4b - Compile FDS MPI release =
#  ======================================

compile_fds_mpi()
{
   # Clean and compile FDS MPI
   cd $FDS_SVNROOT/FDS_Compilation/mpi_intel_linux_64$IB
   make -f ../makefile clean &> /dev/null
   ./make_fds.sh &> $FIREBOT_DIR/output/stage4b
}

check_compile_fds_mpi()
{
   # Check for errors in FDS MPI compilation
   cd $FDS_SVNROOT/FDS_Compilation/mpi_intel_linux_64$IB
   if [ -e "fds_mpi_intel_linux_64$IB" ]
   then
      stage4b_success=true
   else
      echo "Errors from Stage 4b - Compile FDS MPI release:" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage4b >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   # 'performing multi-file optimizations' and 'generating object file' are part of a normal compile
   # grep -v 'feupdateenv ...' ignores a known FDS MPI compiler warning (http://software.intel.com/en-us/forums/showthread.php?t=62806)
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage4b | grep -v 'feupdateenv is not implemented' | grep -v 'performing multi-file optimizations' | grep -v 'generating object file'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 4b - Compile FDS MPI release:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage4b | grep -v 'feupdateenv is not implemented' | grep -v 'performing multi-file optimizations' | grep -v 'generating object file' >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

#  =========================================
#  = Stage 4c - Compile FDS OpenMP release =
#  =========================================

compile_fds_openmp()
{
   # Clean and compile FDS OpenMP
   cd $FDS_SVNROOT/FDS_Compilation/openmp_intel_linux_64
   make -f ../makefile clean &> /dev/null
   ./make_fds.sh &> $FIREBOT_DIR/output/stage4c
}

check_compile_fds_openmp()
{
   # Check for errors in FDS OpenMP compilation
   cd $FDS_SVNROOT/FDS_Compilation/openmp_intel_linux_64
   if [ -e "fds_openmp_intel_linux_64" ]
   then
      stage4c_success=true
   else
      echo "Errors from Stage 4c - Compile FDS OpenMP release:" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage4c >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   # 'performing multi-file optimizations' and 'generating object file' are part of a normal compile
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage4c | grep -v 'performing multi-file optimizations' | grep -v 'generating object file'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 4c - Compile FDS OpenMP release:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage4c | grep -v 'performing multi-file optimizations' | grep -v 'generating object file' >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

#  ======================================
#  = Stage 5pre - Compile SMV utilities =
#  ======================================

compile_smv_utilities()
{  
   # smokeview libraries
   cd $FDS_SVNROOT/SMV/Build/LIBS/lib_linux_intel_64
   echo 'Building Smokeview libraries:' >> $FIREBOT_DIR/output/stage5pre 2>&1
   ./makelibs.sh >> $FIREBOT_DIR/output/stage5pre 2>&1
   echo "" >> $FIREBOT_DIR/output/stage5pre 2>&1

   # smokezip:
   cd $FDS_SVNROOT/Utilities/smokezip/intel_linux_64
   echo 'Compiling smokezip:' >> $FIREBOT_DIR/output/stage5pre 2>&1
   ./make_zip.sh >> $FIREBOT_DIR/output/stage5pre 2>&1
   echo "" >> $FIREBOT_DIR/output/stage5pre 2>&1
   
   # smokediff:
   cd $FDS_SVNROOT/Utilities/smokediff/intel_linux_64
   echo 'Compiling smokediff:' >> $FIREBOT_DIR/output/stage5pre 2>&1
   ./make_diff.sh >> $FIREBOT_DIR/output/stage5pre 2>&1
   echo "" >> $FIREBOT_DIR/output/stage5pre 2>&1
   
   # background:
   cd $FDS_SVNROOT/Utilities/background/intel_linux_32
   echo 'Compiling background:' >> $FIREBOT_DIR/output/stage5pre 2>&1
   ./make_background.sh >> $FIREBOT_DIR/output/stage5pre 2>&1
   echo "" >> $FIREBOT_DIR/output/stage5pre 2>&1

   # wind2fds:
   cd $FDS_SVNROOT/Utilities/wind2fds/intel_linux_64
   echo 'Compiling wind2fds:' >> $FIREBOT_DIR/output/stage5pre 2>&1
   ./make_wind.sh >> $FIREBOT_DIR/output/stage5pre 2>&1
   echo "" >> $FIREBOT_DIR/output/stage5pre 2>&1
}

check_smv_utilities()
{
   # Check for errors in SMV utilities compilation
   cd $FDS_SVNROOT
   if [ -e "$FDS_SVNROOT/Utilities/smokezip/intel_linux_64/smokezip_linux_64" ]  && \
      [ -e "$FDS_SVNROOT/Utilities/smokediff/intel_linux_64/smokediff_linux_64" ]  && \
      [ -e "$FDS_SVNROOT/Utilities/wind2fds/intel_linux_64/wind2fds_linux_64" ]  && \
      [ -e "$FDS_SVNROOT/Utilities/background/intel_linux_32/background" ]
   then
      stage5pre_success=true
   else
      echo "Errors from Stage 5pre - Compile SMV utilities:" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage5pre >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi
}

#  ===================================================
#  = Stage 5 - Run verification cases (release mode) =
#  ===================================================

wait_cases_release_end()
{
   # Scans qstat and waits for verification cases to end
   while [[ `qstat | grep $(whoami)` != '' ]]; do
      JOBS_REMAINING=`qstat | grep $(whoami) | wc -l`
      echo "Waiting for ${JOBS_REMAINING} ${1} cases to complete." >> $FIREBOT_DIR/output/stage5
      TIME_LIMIT_STAGE="5"
      check_time_limit
      sleep 60
   done
}

run_verification_cases_release()
{
   # Start running all FDS verification cases (run all cases on firebot queue)
   cd $FDS_SVNROOT/Verification
   echo 'Running FDS verification cases:' >> $FIREBOT_DIR/output/stage5
   ./Run_FDS_Cases.sh -q $QUEUE >> $FIREBOT_DIR/output/stage5 2>&1
   echo "" >> $FIREBOT_DIR/output/stage5 2>&1

   # Start running all SMV verification cases (run all cases on firebot queue)
   cd $FDS_SVNROOT/Verification/scripts
   echo 'Running SMV verification cases:' >> $FIREBOT_DIR/output/stage5 2>&1
   ./Run_SMV_Cases.sh -q $QUEUE >> $FIREBOT_DIR/output/stage5 2>&1

   # Wait for all verification cases to end
   wait_cases_release_end 'verification'
}

run_validation_cases_release()
{
   #  ===================================
   #  = Run selected FDS validation set =
   #  ===================================

   cd $FDS_SVNROOT/Validation/"$CURRENT_VALIDATION_SET"

   mkdir Current_Results

   # Start running FDS verification cases
   echo "Running FDS validation cases:" >> $FIREBOT_DIR/output/stage5
   echo "Validation Set: ${CURRENT_VALIDATION_SET}" >> $FIREBOT_DIR/output/stage5
   echo "" >> $FIREBOT_DIR/output/stage5 2>&1
   ./Run_All.sh >> $FIREBOT_DIR/output/stage5 2>&1
   echo "" >> $FIREBOT_DIR/output/stage5 2>&1

   # Wait for validation cases to end
   wait_cases_debug_end 'validation'
}

check_cases_release()
{
   # Scan and report any errors in FDS cases
   cd $1

   if [[ `grep -rI 'Run aborted' ${FIREBOT_DIR}/output/stage5` == "" ]] && \
      [[ `grep -rI Segmentation *` == "" ]] && \
      [[ `grep -rI ERROR: *` == "" ]] && \
      [[ `grep -rI 'STOP: Numerical' *` == "" ]] && \
      [[ `grep -rI -A 20 forrtl *` == "" ]]
   then
      stage5_success=true
   else
      grep -rI 'Run aborted' $FIREBOT_DIR/output/stage5 >> $FIREBOT_DIR/output/stage5_errors
      grep -rI Segmentation * >> $FIREBOT_DIR/output/stage5_errors
      grep -rI ERROR: * >> $FIREBOT_DIR/output/stage5_errors
      grep -rI 'STOP: Numerical' * >> $FIREBOT_DIR/output/stage5_errors
      grep -rI -A 20 forrtl * >> $FIREBOT_DIR/output/stage5_errors
      
      echo "Errors from Stage 5 - Run ${2} cases (release mode):" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage5_errors >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi
}

commit_validation_results()
{
   # Copy new FDS files from Current_Results to FDS_Output_Files using Process_Output.csh script for the validation set
   cd $FDS_SVNROOT/Validation/"$CURRENT_VALIDATION_SET"/FDS_Output_Files
   ./Process_Output.csh

   # cd to SVN root
   cd $FDS_SVNROOT

   # Commit new validation results
   svn commit -m "Validationbot: Updated validation results for ${CURRENT_VALIDATION_SET}" &> /dev/null
}

#  ================================
#  = Stage 6a - Compile SMV debug =
#  ================================

compile_smv_db()
{
   # Clean and compile SMV debug
   cd $FDS_SVNROOT/SMV/Build/intel_linux_64_db
   ./make_smv.sh &> $FIREBOT_DIR/output/stage6a
}

check_compile_smv_db()
{
   # Check for errors in SMV debug compilation
   cd $FDS_SVNROOT/SMV/Build/intel_linux_64_db
   if [ -e "smokeview_linux_64_db" ]
   then
      stage6a_success=true
   else
      echo "Errors from Stage 6a - Compile SMV debug:" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage6a >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   # grep -v 'feupdateenv ...' ignores a known FDS MPI compiler warning (http://software.intel.com/en-us/forums/showthread.php?t=62806)
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage6a | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 6a - Compile SMV debug:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage6a | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked' >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

#  =============================================
#  = Stage 6b - Make SMV pictures (debug mode) =
#  =============================================

make_smv_pictures_db()
{
   # Run Make SMV Pictures script (debug mode)
   cd $FDS_SVNROOT/Verification/scripts
   ./Make_SMV_Pictures.sh -d 2>&1 | grep -v FreeFontPath &> $FIREBOT_DIR/output/stage6b
}

check_smv_pictures_db()
{
   # Scan and report any errors in make SMV pictures process
   cd $FIREBOT_DIR
   if [[ `grep -I -E "Segmentation|Error" $FIREBOT_DIR/output/stage6b` == "" ]]
   then
      stage6b_success=true
   else
      cp $FIREBOT_DIR/output/stage6b $FIREBOT_DIR/output/stage6b_errors

      echo "Errors from Stage 6b - Make SMV pictures (debug mode):" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage6b_errors >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi
}

#  ==================================
#  = Stage 6c - Compile SMV release =
#  ==================================

compile_smv()
{
   # Clean and compile SMV
   cd $FDS_SVNROOT/SMV/Build/intel_linux_64
   ./make_smv.sh &> $FIREBOT_DIR/output/stage6c
}

check_compile_smv()
{
   # Check for errors in SMV release compilation
   cd $FDS_SVNROOT/SMV/Build/intel_linux_64
   if [ -e "smokeview_linux_64" ]
   then
      stage6c_success=true
   else
      echo "Errors from Stage 6c - Compile SMV release:" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage6c >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   # grep -v 'feupdateenv ...' ignores a known FDS MPI compiler warning (http://software.intel.com/en-us/forums/showthread.php?t=62806)
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage6c | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 6c - Compile SMV release:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage6c | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked' >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

#  ===============================================
#  = Stage 6d - Make SMV pictures (release mode) =
#  ===============================================

make_smv_pictures()
{
   # Run Make SMV Pictures script (release mode)
   cd $FDS_SVNROOT/Verification/scripts
   ./Make_SMV_Pictures.sh 2>&1 | grep -v FreeFontPath &> $FIREBOT_DIR/output/stage6d
}

check_smv_pictures()
{
   # Scan and report any errors in make SMV pictures process
   cd $FIREBOT_DIR
   if [[ `grep -I -E "Segmentation|Error" $FIREBOT_DIR/output/stage6d` == "" ]]
   then
      stage6d_success=true
   else
      cp $FIREBOT_DIR/output/stage6d $FIREBOT_DIR/output/stage6d_errors

      echo "Errors from Stage 6d - Make SMV pictures (release mode):" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage6d >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi
}

#  ================================
#  = Stage 6e - Make FDS pictures =
#  ================================

make_fds_pictures()
{
   # Run Make FDS Pictures script
   cd $FDS_SVNROOT/Verification
   ./Make_FDS_Pictures.sh &> $FIREBOT_DIR/output/stage6e
}

check_fds_pictures()
{
   # Scan and report any errors in make FDS pictures process
   cd $FIREBOT_DIR
   if [[ `grep -I -E "Segmentation|Error" $FIREBOT_DIR/output/stage6e` == "" ]]
   then
      stage6e_success=true
   else
      cp $FIREBOT_DIR/output/stage6e $FIREBOT_DIR/output/stage6e_errors
      
      echo "Errors from Stage 6e - Make FDS pictures:" >> $ERROR_LOG
      cat $FIREBOT_DIR/output/stage6e_errors >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi
}

#  ====================
#  = Stage 7 - Matlab =
#  ====================

# Functions to check for an available Matlab license

run_matlab_license_test()
{
   # Run simple test to see if Matlab license is available
   cd $FDS_SVNROOT/Utilities/Matlab
   matlab -r "try, disp('Running Matlab License Check'), catch, disp('License Error'), err = lasterror, err.message, err.stack, end, exit" &> $FIREBOT_DIR/output/stage7_matlab_license
}

scan_matlab_license_test()
{
   # Check for failed license
   if [[ `grep "License checkout failed" $FIREBOT_DIR/output/stage7_matlab_license` == "" ]]
   then
      # Continue along
      :
   else
      TIME_LIMIT_STAGE="7"
      check_time_limit
      # Wait 5 minutes until retry
      sleep 300
      check_matlab_license_server
   fi
}

check_matlab_license_server()
{
   run_matlab_license_test
   scan_matlab_license_test
}

#  ============================================================
#  = Stage 7a - Matlab plotting and statistics (verification) =
#  ============================================================

run_matlab_verification()
{
   # Run Matlab plotting script
   cd $FDS_SVNROOT/Utilities/Matlab
   matlab -r "try, disp('Running Matlab Verification script'), FDS_verification_script, catch, disp('Error'), err = lasterror, err.message, err.stack, end, exit" &> $FIREBOT_DIR/output/stage7a_verification
}

check_matlab_verification()
{
   # Scan and report any errors in Matlab scripts
   cd $FIREBOT_DIR
   if [[ `grep -A 50 "Error" $FIREBOT_DIR/output/stage7a_verification` == "" ]]
   then
      stage7a_success=true
   else
      echo "Warnings from Stage 7a - Matlab plotting and statistics (verification):" >> $WARNING_LOG
      grep -A 50 "Error" $FIREBOT_DIR/output/stage7a_verification >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

check_verification_stats()
{
   # Check for existence of verification statistics output file
   cd $FDS_SVNROOT/Utilities/Matlab
   if [ -e "FDS_verification_scatterplot_output.csv" ]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 7a - Matlab plotting and statistics (verification):" >> $WARNING_LOG
      echo "Error: The verification statistics output file does not exist." >> $WARNING_LOG
      echo "Expected the file Utilities/Matlab/FDS_verification_scatterplot_output.csv" >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi

   # Scan and report warnings for any verification cases that are outside of their specified error tolerance
   cd $FDS_SVNROOT/Utilities/Matlab
   if [[ `grep "Out of Tolerance" FDS_verification_scatterplot_output.csv` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 7a - Matlab plotting and statistics (verification):" >> $WARNING_LOG
      echo "The following cases are outside of their specified error tolerance:" >> $WARNING_LOG
      grep "Out of Tolerance" FDS_verification_scatterplot_output.csv >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi

   # Scan and report any case warnings in Matlab scripts
   cd $FIREBOT_DIR
   if [[ `grep "Matlab Warning" $FIREBOT_DIR/output/stage7a_verification` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 7a - Matlab plotting and statistics (verification):" >> $WARNING_LOG
      grep "Matlab Warning" $FIREBOT_DIR/output/stage7a_verification >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

#  ===========================================
#  = Stage 7b - Matlab plotting (validation) =
#  ===========================================

run_matlab_validation()
{
   # Run Matlab plotting script
   cd $FDS_SVNROOT/Utilities/Matlab
   matlab -r "try, disp('Running Matlab Validation script'), FDS_validation_script, catch, disp('Error'), err = lasterror, err.message, err.stack, end, exit" &> $FIREBOT_DIR/output/stage7b_validation
}

check_matlab_validation()
{
   # Scan and report any errors in Matlab scripts
   cd $FIREBOT_DIR
   if [[ `grep -A 50 "Error" $FIREBOT_DIR/output/stage7b_validation` == "" ]]
   then
      stage7b_success=true
   else
      echo "Warnings from Stage 7b - Matlab plotting (validation):" >> $WARNING_LOG
      grep -A 50 "Error" $FIREBOT_DIR/output/stage7b_validation >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

#  ======================================
#  = Stage 7c - FDS run time statistics =
#  ======================================

generate_timing_stats()
{
   cd $FDS_SVNROOT/Utilities/Scripts
   ./fds_timing_stats.sh
}

archive_timing_stats()
{
   cd $FDS_SVNROOT/Utilities/Scripts
   cp fds_timing_stats.csv "$FIREBOT_DIR/history/${SVN_REVISION}_timing.csv"
}

#  ==================================
#  = Stage 8 - Build FDS-SMV guides =
#  ==================================

check_guide()
{
   # Scan and report any errors or warnings in build process for guides
   cd $FIREBOT_DIR
   if [[ `grep -I "successfully" $1` == "" ]]
   then
      # There were errors/warnings in the guide build process
      echo "Warnings from Stage 8 - Build FDS-SMV Guides:" >> $WARNING_LOG
      echo $3 >> $WARNING_LOG # Name of guide
      cat $1 >> $WARNING_LOG # Contents of log file
      echo "" >> $WARNING_LOG
   else
      # Guide built successfully; there were no errors/warnings
      # Copy guide to Firebot's local website
      cp $2 /var/www/html/firebot/manuals/
   fi
}

make_fds_user_guide()
{
   cd $FDS_SVNROOT/Manuals/FDS_User_Guide

   # Build FDS User Guide
   ./make_guide.sh &> $FIREBOT_DIR/output/stage8_fds_user_guide

   # Check guide for completion and copy to website if successful
   check_guide $FIREBOT_DIR/output/stage8_fds_user_guide $FDS_SVNROOT/Manuals/FDS_User_Guide/FDS_User_Guide.pdf 'FDS User Guide'
}

make_fds_technical_guide()
{
   cd $FDS_SVNROOT/Manuals/FDS_Technical_Reference_Guide

   # Build FDS Technical Guide
   ./make_guide.sh &> $FIREBOT_DIR/output/stage8_fds_technical_guide

   # Check guide for completion and copy to website if successful
   check_guide $FIREBOT_DIR/output/stage8_fds_technical_guide $FDS_SVNROOT/Manuals/FDS_Technical_Reference_Guide/FDS_Technical_Reference_Guide.pdf 'FDS Technical Reference Guide'
}

make_fds_verification_guide()
{
   cd $FDS_SVNROOT/Manuals/FDS_Verification_Guide

   # Build FDS Verification Guide
   ./make_guide.sh &> $FIREBOT_DIR/output/stage8_fds_verification_guide

   # Check guide for completion and copy to website if successful
   check_guide $FIREBOT_DIR/output/stage8_fds_verification_guide $FDS_SVNROOT/Manuals/FDS_Verification_Guide/FDS_Verification_Guide.pdf 'FDS Verification Guide'
}

make_fds_validation_guide()
{
   cd $FDS_SVNROOT/Manuals/FDS_Validation_Guide

   # Build FDS Validation Guide
   ./make_guide.sh &> $FIREBOT_DIR/output/stage8_fds_validation_guide

   # Check guide for completion and copy to website if successful
   check_guide $FIREBOT_DIR/output/stage8_fds_validation_guide $FDS_SVNROOT/Manuals/FDS_Validation_Guide/FDS_Validation_Guide.pdf 'FDS Validation Guide'
}

make_fds_configuration_management_plan()
{
   cd $FDS_SVNROOT/Manuals/FDS_Configuration_Management_Plan

   # Build FDS Configuration Management Plan
   ./make_guide.sh &> $FIREBOT_DIR/output/stage8_fds_configuration_management_plan

   # Check guide for completion and copy to website if successful
   check_guide $FIREBOT_DIR/output/stage8_fds_configuration_management_plan $FDS_SVNROOT/Manuals/FDS_Configuration_Management_Plan/FDS_Configuration_Management_Plan.pdf 'FDS Configuration Management Plan'
}

make_smv_user_guide()
{
   cd $FDS_SVNROOT/Manuals/SMV_User_Guide

   # Build SMV User Guide
   ./make_guide.sh &> $FIREBOT_DIR/output/stage8_smv_user_guide

   # Check guide for completion and copy to website if successful
   check_guide $FIREBOT_DIR/output/stage8_smv_user_guide $FDS_SVNROOT/Manuals/SMV_User_Guide/SMV_User_Guide.pdf 'SMV User Guide'
}

make_smv_technical_guide()
{
   cd $FDS_SVNROOT/Manuals/SMV_Technical_Reference_Guide

   # Build SMV Technical Guide
   ./make_guide.sh &> $FIREBOT_DIR/output/stage8_smv_technical_guide

   # Check guide for completion and copy to website if successful
   check_guide $FIREBOT_DIR/output/stage8_smv_technical_guide $FDS_SVNROOT/Manuals/SMV_Technical_Reference_Guide/SMV_Technical_Reference_Guide.pdf 'SMV Technical Reference Guide'
}

make_smv_verification_guide()
{
   cd $FDS_SVNROOT/Manuals/SMV_Verification_Guide

   # Build SMV Verification Guide
   ./make_guide.sh &> $FIREBOT_DIR/output/stage8_smv_verification_guide

   # Check guide for completion and copy to website if successful
   check_guide $FIREBOT_DIR/output/stage8_smv_verification_guide $FDS_SVNROOT/Manuals/SMV_Verification_Guide/SMV_Verification_Guide.pdf 'SMV Verification Guide'
}

#  =====================================================
#  = Build status reporting - email and save functions =
#  =====================================================

save_build_status()
{
   cd $FIREBOT_DIR
   # Save status outcome of build to a text file
   if [[ -e $WARNING_LOG && -e $ERROR_LOG ]]
   then
     echo "" >> $ERROR_LOG
     cat $WARNING_LOG >> $ERROR_LOG
     echo "Build failure and warnings for Revision ${SVN_REVISION}." > "$FIREBOT_DIR/history/${SVN_REVISION}.txt"
     cat $ERROR_LOG > "$FIREBOT_DIR/history/${SVN_REVISION}_errors.txt"

   # Check for errors only
   elif [ -e $ERROR_LOG ]
   then
      echo "Build failure for Revision ${SVN_REVISION}." > "$FIREBOT_DIR/history/${SVN_REVISION}.txt"
      cat $ERROR_LOG > "$FIREBOT_DIR/history/${SVN_REVISION}_errors.txt"

   # Check for warnings only
   elif [ -e $WARNING_LOG ]
   then
      echo "Revision ${SVN_REVISION} has warnings." > "$FIREBOT_DIR/history/${SVN_REVISION}.txt"
      cat $WARNING_LOG > "$FIREBOT_DIR/history/${SVN_REVISION}_warnings.txt"

   # No errors or warnings
   else
      echo "Build success! Revision ${SVN_REVISION} passed all build tests." > "$FIREBOT_DIR/history/${SVN_REVISION}.txt"
   fi
}

email_build_status()
{
   cd $FIREBOT_DIR
   # Check for warnings and errors
   if [[ -e $WARNING_LOG && -e $ERROR_LOG ]]
   then
     # Send email with failure message and warnings, body of email contains appropriate log file
     mail -s "[${1}@$hostname] ${2} failure and warnings for Revision ${SVN_REVISION}." $mailToFDS < $ERROR_LOG > /dev/null

   # Check for errors only
   elif [ -e $ERROR_LOG ]
   then
      # Send email with failure message, body of email contains error log file
      mail -s "[${1}@$hostname] ${2} failure for Revision ${SVN_REVISION}." $mailToFDS < $ERROR_LOG > /dev/null

   # Check for warnings only
   elif [ -e $WARNING_LOG ]
   then
      # Send email with success message, include warnings
      mail -s "[${1}@$hostname] ${2} success, with warnings. Revision ${SVN_REVISION} passed all build tests." $mailToFDS < $WARNING_LOG > /dev/null

   # No errors or warnings
   else
      # Send success message with links to nightly manuals
      stop_time=`date`
      echo "-------------------------------" >> $TIME_LOG
      echo "Host OS: Linux " >> $TIME_LOG
      echo "Host Name: $hostname " >> $TIME_LOG
      if [ $FIREBOT_MODE == "validation" ] ; then
         echo "Validation Set: ${CURRENT_VALIDATION_SET} " >> $TIME_LOG
      fi
      echo "Start Time: $start_time " >> $TIME_LOG
      echo "Stop Time: $stop_time " >> $TIME_LOG
      echo "-------------------------------" >> $TIME_LOG
      echo "Nightly Manuals (private):  http://blaze.nist.gov/firebot" >> $TIME_LOG
      echo "Nightly Manuals (public):   https://drive.google.com/folderview?id=0B_wB1pJL2bFQaDJaOFNnUDR4LXM#list" >> $TIME_LOG
      echo "-------------------------------" >> $TIME_LOG
      mail -s "[${1}@$hostname] ${2} success! Revision ${SVN_REVISION} passed all build tests." $mailToFDS < $TIME_LOG > /dev/null
   fi
}

#  ============================
#  = Primary script execution =
#  ============================

hostname=`hostname`
start_time=`date`

### Clean up on start ###
clean_firebot_history

### Stage 0 ###
update_and_compile_cfast

### Stage 1 ###
clean_svn_repo
do_svn_checkout
check_svn_checkout
# Only run if -s option (skip SVN properties) is not used
if [[ ! $SKIP_SVN_PROPS ]] ; then
   fix_svn_properties
fi

### Stage 2a ###
compile_fds_db
check_compile_fds_db

### Stage 2b ###
compile_fds_mpi_db
check_compile_fds_mpi_db

### Stage 2c ###
compile_fds_openmp_db
check_compile_fds_openmp_db

# Depends on successful FDS OpenMP debug compile
# if [[ $stage2c_success ]] ; then
   # inspect_fds_openmp_db
   # check_inspect_fds_openmp_db
# fi

### Stage 3 ###
# Only run if firebot is in "validation" mode
if [ $FIREBOT_MODE == "validation" ] ; then
   select_validation_set
fi

# Depends on successful FDS debug compile
if [[ $stage2a_success && $stage2b_success && $FIREBOT_MODE == "verification" ]] ; then
   run_verification_cases_debug
   check_cases_debug $FDS_SVNROOT/Verification 'verification'
elif [[ $stage2a_success && $stage2b_success && $FIREBOT_MODE == "validation" ]] ; then
   run_validation_cases_debug
   check_cases_debug $FDS_SVNROOT/Validation 'validation'
fi

### Stage 4a ###
compile_fds
check_compile_fds

### Stage 4b ###
compile_fds_mpi
check_compile_fds_mpi

### Stage 4c ###
compile_fds_openmp
check_compile_fds_openmp

### Stage 5pre ###
# Only run if firebot is in "verification" mode
if [ $FIREBOT_MODE == "verification" ] ; then
   compile_smv_utilities
   check_smv_utilities
fi

### Stage 5 ###
# Depends on successful FDS compile
if [[ $stage4a_success && $stage4b_success && $FIREBOT_MODE == "verification" ]] ; then
   run_verification_cases_release
   check_cases_release $FDS_SVNROOT/Verification 'verification'
elif [[ $stage4a_success && $stage4b_success && $FIREBOT_MODE == "validation" ]] ; then
   run_validation_cases_release
   check_cases_release $FDS_SVNROOT/Validation 'validation'
fi

# Depends on successful run of validation cases
if [[ $stage5_success && $FIREBOT_MODE == "validation" ]] ; then
   commit_validation_results
fi

#  ======================================================================
#  = Only run the following stages if firebot is in "verification" mode =
#  ======================================================================

if [ $FIREBOT_MODE == "verification" ] ; then
   ### Stage 6a ###
   compile_smv_db
   check_compile_smv_db

   ### Stage 6b ###
   # Depends on successful SMV debug compile
   if [[ $stage6a_success ]] ; then
      make_smv_pictures_db
      check_smv_pictures_db
   fi

   ### Stage 6c ###
   compile_smv
   check_compile_smv

   ### Stage 6d ###
   # Depends on successful SMV compile
   if [[ $stage6c_success ]] ; then
      make_smv_pictures
      check_smv_pictures
   fi

   ### Stage 6e ###
   # Depends on successful SMV compile
   if [[ $stage6c_success ]] ; then
      make_fds_pictures
      check_fds_pictures
   fi

   ### Stage 7a ###
   check_matlab_license_server
   run_matlab_verification
   check_matlab_verification
   check_verification_stats

   ### Stage 7b ###
   check_matlab_license_server
   run_matlab_validation
   check_matlab_validation

   ### Stage 7c ###
   generate_timing_stats
   archive_timing_stats

   ### Stage 8 ###
   make_fds_user_guide
   make_fds_verification_guide
   make_fds_technical_guide
   make_fds_validation_guide
   make_smv_user_guide
   make_smv_technical_guide
   make_smv_verification_guide
   make_fds_configuration_management_plan
fi

### Wrap up and report results ###
set_files_world_readable
if [ $FIREBOT_MODE == "verification" ] ; then
   save_build_status
   email_build_status 'Firebot' 'Build'
elif [ $FIREBOT_MODE == "validation" ] ; then
   email_build_status 'Validationbot' 'Validation'
fi

