#!/bin/bash

# Firebot
# FDS automatIc veRification and validation tEst bot
# Kristopher Overholt
# 6/22/2012

#  ===================
#  = Input variables =
#  ===================

mailTo="kevin.mcgrattan@nist.gov, randall.mcdermott@nist.gov, glenn.forney@nist.gov, craig.weinschenk@nist.gov, kristopher.overholt@nist.gov"
SVNROOT="/home/firebot/FDS-SMV"
FIREBOT_DIR="/home/firebot/firebot"
SVN_REVISION=$1

#  =========================
#  = External dependencies =
#  =========================
#
#   This script expects the following dependencies to be in place:
#   
#   cfast (for Stage 5 - Run_SMV_Cases.sh):
#      ~/cfast/CFAST/intel_linux_64/cfast6_linux_64
#

#  ====================
#  = End user warning =
#  ====================

# Warn if running as user other than firebot
if [[ `whoami` == "firebot" ]];
   then
      # Continue along
      :
   else
      echo "Warning: You are running the Firebot script as an end user."
      echo "This script can modify and erase your repository."
      echo "If you wish to continue, edit the script and remove this warning."
      echo "Terminating script."
      exit
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
         echo 'Sending email'
         echo -e "Firebot has been running for more than 12 hours in Stage ${TIME_LIMIT_STAGE}. \n\nPlease ensure that there are no problems. \n\nThis is a notification only and does not terminate Firebot." | mail -s "[Firebot] Notice: Firebot has been running for more than 12 hours." $mailTo > /dev/null
         TIME_LIMIT_EMAIL_NOTIFICATION="sent"
      fi
   fi
}

#  ========================
#  ========================
#  = Firebot Build Stages =
#  ========================
#  ========================

#  ============================
#  = Stage 1 - SVN operations =
#  ============================

clean_svn_repo()
{
   # Initialize and start with fresh repo
   # Clean Firebot metafiles
   cd $FIREBOT_DIR
   rm output/*

   # Clean up temporary unversioned and modified versioned repository files
   cd $SVNROOT
   svn revert -Rq *
   svn status --no-ignore | grep '^\?' | sed 's/^\?      //'  | xargs -Ixx rm -rf xx
}

do_svn_checkout()
{
   # If an SVN revision number is specified, then get that revision
   if [[ $SVN_REVISION != "" ]]; then
      echo "Checking out revision r${SVN_REVISION}." > $FIREBOT_DIR/output/stage1
      svn update -r $SVN_REVISION >> $FIREBOT_DIR/output/stage1 2>&1
   # If no SVN revision number is specified, then get the latest revision
   else
      echo "Checking out latest revision." > $FIREBOT_DIR/output/stage1
      svn update >> $FIREBOT_DIR/output/stage1 2>&1
      SVN_REVISION=`tail -n 1 $FIREBOT_DIR/output/stage1 | sed "s/[^0-9]//g"`
   fi
}

check_svn_checkout()
{
   cd $SVNROOT
   # Check for SVN errors
   if [[ `grep -E 'Updated|At revision' $FIREBOT_DIR/output/stage1 | wc -l` -ne 1 ]];
   then
      BUILD_STAGE_FAILURE="Stage 1: SVN Operations"
      ERROR_LOG=$FIREBOT_DIR/output/stage1
      save_build_status
      email_error_message
   else
      # Continue along
      :
   fi
}

#  =============================
#  = Stage 2a - Compile FDS DB =
#  =============================

compile_fds_db()
{
   # Clean and compile FDS DB
   cd $SVNROOT/FDS_Compilation/intel_linux_64_db
   make --makefile ../makefile clean &> /dev/null
   ./make_fds.sh &> $FIREBOT_DIR/output/stage2a
}

check_compile_fds_db()
{
   # Check for errors in FDS DB compilation
   cd $SVNROOT/FDS_Compilation/intel_linux_64_db
   if [ -e "fds_intel_linux_64_db" ]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 2a: FDS DB Compilation"
      ERROR_LOG=$FIREBOT_DIR/output/stage2a
      save_build_status
      email_error_message
   fi

   # Check for compiler warnings/remarks
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2a` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 2a warnings:" >> $FIREBOT_DIR/output/warnings
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2a >> $FIREBOT_DIR/output/warnings
      echo "" >> $FIREBOT_DIR/output/warnings
   fi
}

#  =================================
#  = Stage 2b - Compile FDS MPI DB =
#  =================================

compile_fds_mpi_db()
{
   # Clean and compile FDS MPI DB
   cd $SVNROOT/FDS_Compilation/mpi_intel_linux_64_db
   make --makefile ../makefile clean &> /dev/null
   ./make_fds.sh &> $FIREBOT_DIR/output/stage2b
}

check_compile_fds_mpi_db()
{
   # Check for errors in FDS MPI DB compilation
   cd $SVNROOT/FDS_Compilation/mpi_intel_linux_64_db
   if [ -e "fds_mpi_intel_linux_64_db" ]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 2b: FDS MPI DB Compilation"
      ERROR_LOG=$FIREBOT_DIR/output/stage2b
      save_build_status
      email_error_message
   fi

   # Check for compiler warnings/remarks
   # grep -v 'feupdateenv ...' ignores a known FDS MPI compiler warning (http://software.intel.com/en-us/forums/showthread.php?t=62806)
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2b | grep -v 'feupdateenv is not implemented'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 2b warnings:" >> $FIREBOT_DIR/output/warnings
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2b | grep -v 'feupdateenv is not implemented' >> $FIREBOT_DIR/output/warnings
      echo "" >> $FIREBOT_DIR/output/warnings
   fi
}

#  ================================================
#  = Stage 3 - Run verification cases (short run) =
#  ================================================

wait_verification_cases_short_start()
{
   # Scans qstat and waits for verification cases to start
   while [[ `qstat | grep $(whoami) | grep Q` != '' ]]; do
      JOBS_REMAINING=`qstat | grep $(whoami) | grep Q | wc -l`
      echo "Waiting for ${JOBS_REMAINING} verification cases to start." >> $FIREBOT_DIR/output/stage3
      TIME_LIMIT_STAGE="3"
      check_time_limit
      sleep 30
   done
}

wait_verification_cases_short_end()
{
   # Scans qstat and waits for verification cases to end
   while [[ `qstat | grep $(whoami)` != '' ]]; do
      JOBS_REMAINING=`qstat | grep $(whoami) | wc -l`
      echo "Waiting for ${JOBS_REMAINING} verification cases to complete." >> $FIREBOT_DIR/output/stage3
      TIME_LIMIT_STAGE="3"
      check_time_limit
      sleep 30
   done
}

run_verification_cases_short()
{

   #  ============================
   #  = Run all FDS serial cases =
   #  ============================

   cd $SVNROOT/Verification

   # Submit FDS verification cases and wait for them to start (run serial cases in debug mode on fire70s queue)
   ./Run_FDS_Cases.sh -c serial -d -q fire70s &> $FIREBOT_DIR/output/stage3
   wait_verification_cases_short_start

   # Wait some additional time for all cases to start
   sleep 30

   # Stop all cases
   ./Run_FDS_Cases.sh -c serial -d -s >> $FIREBOT_DIR/output/stage3 2>&1

   # Wait for serial verification cases to end
   wait_verification_cases_short_end

   #  =========================
   #  = Run all FDS MPI cases =
   #  =========================

   cd $SVNROOT/Verification

   # Submit FDS verification cases and wait for them to start (run MPI cases in debug mode on fire70s queue)
   ./Run_FDS_Cases.sh -c mpi -d -q fire70s &> $FIREBOT_DIR/output/stage3
   wait_verification_cases_short_start

   # Wait some additional time for all cases to start
   sleep 30

   # Stop all cases
   ./Run_FDS_Cases.sh -c mpi -d -s >> $FIREBOT_DIR/output/stage3 2>&1

   # Wait for serial verification cases to end
   wait_verification_cases_short_end

   #  =====================
   #  = Run all SMV cases =
   #  =====================

   # cd $SVNROOT/Verification/scripts

   # # Submit SMV verification cases and wait for them to start (run all cases in debug mode on fire70s queue)
   # ./scripts/Run_SMV_Cases.sh -d -q fire70s >> $FIREBOT_DIR/output/stage3 2>&1
   # wait_verification_cases_short_start

   # # Wait some additional time for all cases to start
   # sleep 30

   # # Stop all cases
   # ./scripts/Run_SMV_Cases.sh -d -s >> $FIREBOT_DIR/output/stage3 2>&1

   # # Wait for SMV verification cases to end
   # wait_verification_cases_short_end

   #  ======================
   #  = Remove .stop files =
   #  ======================

   # Remove all .stop files from Verification directories (recursively)
   cd $SVNROOT/Verification
   find . -name '*.stop' -exec rm -f {} \;
}

check_verification_cases_short()
{
   # Scan and report any errors in FDS verification cases
   cd $SVNROOT/Verification

   if [[ `grep 'Run aborted' -rI ${FIREBOT_DIR}/output/stage3` == "" ]] && \
      [[ `grep ERROR: -rI *` == "" ]] && \
      [[ `grep 'STOP: Numerical' -rI *` == "" ]] && \
      [[ `grep -A 20 forrtl -rI *` == "" ]]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 3: FDS Verification Cases"
      
      grep 'Run aborted' -rI $FIREBOT_DIR/output/stage3 >> $FIREBOT_DIR/output/stage3_errors
      grep ERROR: -rI * >> $FIREBOT_DIR/output/stage3_errors
      grep 'STOP: Numerical' -rI * >> $FIREBOT_DIR/output/stage3_errors
      grep -A 20 forrtl -rI * >> $FIREBOT_DIR/output/stage3_errors
      
      ERROR_LOG=$FIREBOT_DIR/output/stage3_errors
      save_build_status
      email_error_message
   fi
}

#  ==================================
#  = Stage 4a - Compile FDS release =
#  ==================================

compile_fds()
{
   # Clean and compile FDS
   cd $SVNROOT/FDS_Compilation/intel_linux_64
   make --makefile ../makefile clean &> /dev/null
   ./make_fds.sh &> $FIREBOT_DIR/output/stage4a
}

check_compile_fds()
{
   # Check for errors in FDS compilation
   cd $SVNROOT/FDS_Compilation/intel_linux_64
   if [ -e "fds_intel_linux_64" ]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 4a: FDS Release Compilation"
      ERROR_LOG=$FIREBOT_DIR/output/stage4a
      save_build_status
      email_error_message
   fi

   # Check for compiler warnings/remarks
   # 'performing multi-file optimizations' and 'generating object file' are part of a normal compile
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage4a | grep -v 'performing multi-file optimizations' | grep -v 'generating object file'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 4a warnings:" >> $FIREBOT_DIR/output/warnings
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage4a | grep -v 'performing multi-file optimizations' | grep -v 'generating object file' >> $FIREBOT_DIR/output/warnings
      echo "" >> $FIREBOT_DIR/output/warnings
   fi
}

#  ======================================
#  = Stage 4b - Compile FDS MPI release =
#  ======================================

compile_fds_mpi()
{
   # Clean and compile FDS MPI
   cd $SVNROOT/FDS_Compilation/mpi_intel_linux_64
   make --makefile ../makefile clean &> /dev/null
   ./make_fds.sh &> $FIREBOT_DIR/output/stage4b
}

check_compile_fds_mpi()
{
   # Check for errors in FDS MPI compilation
   cd $SVNROOT/FDS_Compilation/mpi_intel_linux_64
   if [ -e "fds_mpi_intel_linux_64" ]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 4b: FDS MPI Release Compilation"
      ERROR_LOG=$FIREBOT_DIR/output/stage4b
      save_build_status
      email_error_message
   fi

   # Check for compiler warnings/remarks
   # 'performing multi-file optimizations' and 'generating object file' are part of a normal compile
   # grep -v 'feupdateenv ...' ignores a known FDS MPI compiler warning (http://software.intel.com/en-us/forums/showthread.php?t=62806)
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage4b | grep -v 'feupdateenv is not implemented' | grep -v 'performing multi-file optimizations' | grep -v 'generating object file'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 4b warnings:" >> $FIREBOT_DIR/output/warnings
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage4b | grep -v 'feupdateenv is not implemented' | grep -v 'performing multi-file optimizations' | grep -v 'generating object file' >> $FIREBOT_DIR/output/warnings
      echo "" >> $FIREBOT_DIR/output/warnings
   fi
}

#  ===============================================
#  = Stage 5 - Run verification cases (long run) =
#  ===============================================

wait_verification_cases_long_end()
{
   # Scans qstat and waits for verification cases to end
   while [[ `qstat | grep $(whoami)` != '' ]]; do
      JOBS_REMAINING=`qstat | grep $(whoami) | wc -l`
      echo "Waiting for ${JOBS_REMAINING} verification cases to complete." >> $FIREBOT_DIR/output/stage5
      TIME_LIMIT_STAGE="5"
      check_time_limit
      sleep 60
   done
}

run_verification_cases_long()
{
   # Start running all FDS verification cases (run all cases on fire70s queue)
   cd $SVNROOT/Verification
   echo 'Running FDS verification cases:' > $FIREBOT_DIR/output/stage5
   ./Run_FDS_Cases.sh -q fire70s >> $FIREBOT_DIR/output/stage5 2>&1
   echo "" >> $FIREBOT_DIR/output/stage5 2>&1

   # Start running all SMV verification cases (run all cases on fire70s queue)
   cd $SVNROOT/Verification/scripts
   echo 'Running SMV verification cases:' >> $FIREBOT_DIR/output/stage5 2>&1
   ./Run_SMV_Cases.sh -q fire70s >> $FIREBOT_DIR/output/stage5 2>&1

   # Wait for all verification cases to end
   wait_verification_cases_long_end
}

check_verification_cases_long()
{
   # Scan and report any errors in FDS verification cases
   cd $SVNROOT/Verification

   if [[ `grep 'Run aborted' -rI ${FIREBOT_DIR}/output/stage5` == "" ]] && \
      [[ `grep ERROR: -rI *` == "" ]] && \
      [[ `grep 'STOP: Numerical' -rI *` == "" ]] && \
      [[ `grep -A 20 forrtl -rI *` == "" ]]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 5: FDS-SMV Verification Cases"
      
      grep 'Run aborted' -rI $FIREBOT_DIR/output/stage5 >> $FIREBOT_DIR/output/stage5_errors
      grep ERROR: -rI * >> $FIREBOT_DIR/output/stage5_errors
      grep 'STOP: Numerical' -rI * >> $FIREBOT_DIR/output/stage5_errors
      grep -A 20 forrtl -rI * >> $FIREBOT_DIR/output/stage5_errors
      
      ERROR_LOG=$FIREBOT_DIR/output/stage5_errors
      save_build_status
      email_error_message
   fi
}

#  ====================================
#  = Stage 6a - Compile SMV utilities =
#  ====================================

compile_smv_utilities()
{  
   # smokezip:
   cd $SVNROOT/Utilities/smokezip/intel_linux_64
   echo 'Compiling smokezip:' > $FIREBOT_DIR/output/stage6a 2>&1
   ./make_zip.sh >> $FIREBOT_DIR/output/stage6a 2>&1
   echo "" >> $FIREBOT_DIR/output/stage6a 2>&1
   
   # smokediff:
   cd $SVNROOT/Utilities/smokediff/intel_linux_64
   echo 'Compiling smokediff:' >> $FIREBOT_DIR/output/stage6a 2>&1
   ./make_diff.sh >> $FIREBOT_DIR/output/stage6a 2>&1
   echo "" >> $FIREBOT_DIR/output/stage6a 2>&1
   
   # background:
   cd $SVNROOT/Utilities/background/intel_linux_32
   echo 'Compiling background:' >> $FIREBOT_DIR/output/stage6a 2>&1
   ./make_background.sh >> $FIREBOT_DIR/output/stage6a 2>&1
}

check_smv_utilities()
{
   # Check for errors in SMV utilities compilation
   cd $SVNROOT
   if [ -e "$SVNROOT/Utilities/smokezip/intel_linux_64/smokezip_linux_64" ]  && \
      [ -e "$SVNROOT/Utilities/smokediff/intel_linux_64/smokediff_linux_64" ]  && \
      [ -e "$SVNROOT/Utilities/background/intel_linux_32/background" ]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 6a: SMV Utilities Compilation"
      ERROR_LOG=$FIREBOT_DIR/output/stage6a
      save_build_status
      email_error_message
   fi
}

#  ==================================
#  = Stage 6b - Compile SMV test DB =
#  ==================================

compile_smv_test_db()
{
   # Clean and compile SMV test DB
   cd $SVNROOT/SMV/Build/intel_linux_test_64_dbg
   make --makefile ../Makefile clean &> /dev/null
   ./make_smv.sh &> $FIREBOT_DIR/output/stage6b
}

check_compile_smv_test_db()
{
   # Check for errors in SMV test DB compilation
   cd $SVNROOT/SMV/Build/intel_linux_test_64_dbg
   if [ -e "smokeview_linux_test_64_dbg" ]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 6b: SMV Test DB Compilation"
      ERROR_LOG=$FIREBOT_DIR/output/stage6b
      save_build_status
      email_error_message
   fi

   # Check for compiler warnings/remarks
   # grep -v 'feupdateenv ...' ignores a known FDS MPI compiler warning (http://software.intel.com/en-us/forums/showthread.php?t=62806)
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage6b | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 6b warnings:" >> $FIREBOT_DIR/output/warnings
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage6b | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked' >> $FIREBOT_DIR/output/warnings
      echo "" >> $FIREBOT_DIR/output/warnings
   fi
}

#  ==================================================
#  = Stage 6c - Make SMV pictures (test debug mode) =
#  ==================================================

make_smv_pictures_db()
{
   # Run Make SMV Pictures script (test debug mode)
   cd $SVNROOT/Verification/scripts
   ./Make_SMV_Pictures.sh -d &> $FIREBOT_DIR/output/stage6c
}

check_smv_pictures_db()
{
   # Scan and report any errors in make SMV pictures process
   cd $FIREBOT_DIR
   if [[ `grep -B 50 -A 50 "Segmentation" -I $FIREBOT_DIR/output/stage6c` == "" ]]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 6c: Make SMV Pictures (Test Debug Mode)"
      grep -B 50 -A 50 "Segmentation" -I $FIREBOT_DIR/output/stage6c > $FIREBOT_DIR/output/stage6c_errors
      ERROR_LOG=$FIREBOT_DIR/output/stage6c_errors
      save_build_status
      email_error_message
   fi
}

#  ===============================
#  = Stage 6d - Compile SMV test =
#  ===============================

compile_smv_test()
{
   # Clean and compile SMV DB
   cd $SVNROOT/SMV/Build/intel_linux_test_64
   make --makefile ../Makefile clean &> /dev/null
   ./make_smv.sh &> $FIREBOT_DIR/output/stage6d
}

check_compile_smv_test()
{
   # Check for errors in SMV test compilation
   cd $SVNROOT/SMV/Build/intel_linux_test_64
   if [ -e "smokeview_linux_test_64" ]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 6d: SMV Test Compilation"
      ERROR_LOG=$FIREBOT_DIR/output/stage6d
      save_build_status
      email_error_message
   fi

   # Check for compiler warnings/remarks
   # grep -v 'feupdateenv ...' ignores a known FDS MPI compiler warning (http://software.intel.com/en-us/forums/showthread.php?t=62806)
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage6d | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 6d warnings:" >> $FIREBOT_DIR/output/warnings
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage6d | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked' >> $FIREBOT_DIR/output/warnings
      echo "" >> $FIREBOT_DIR/output/warnings
   fi
}

#  ============================================
#  = Stage 6e - Make SMV pictures (test mode) =
#  ============================================

make_smv_pictures_test()
{
   # Run Make SMV Pictures script (test mode)
   cd $SVNROOT/Verification/scripts
   ./Make_SMV_Pictures.sh &> $FIREBOT_DIR/output/stage6e
}

check_smv_pictures_test()
{
   # Scan and report any errors in make SMV pictures process
   cd $FIREBOT_DIR
   if [[ `grep -B 50 -A 50 "Segmentation" -I $FIREBOT_DIR/output/stage6e` == "" ]]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 6e: Make SMV Pictures (Test Mode)"
      grep -B 50 -A 50 "Segmentation" -I $FIREBOT_DIR/output/stage6e > $FIREBOT_DIR/output/stage6e_errors
      ERROR_LOG=$FIREBOT_DIR/output/stage6e_errors
      save_build_status
      email_error_message
   fi
}

#  ==================================
#  = Stage 6f - Compile SMV release =
#  ==================================

compile_smv()
{
   # Clean and compile SMV
   cd $SVNROOT/SMV/Build/intel_linux_64
   make --makefile ../Makefile clean &> /dev/null
   ./make_smv.sh &> $FIREBOT_DIR/output/stage6f
}

check_compile_smv()
{
   # Check for errors in SMV release compilation
   cd $SVNROOT/SMV/Build/intel_linux_64
   if [ -e "smokeview_linux_64" ]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 6f: SMV Release Compilation"
      ERROR_LOG=$FIREBOT_DIR/output/stage6f
      save_build_status
      email_error_message
   fi

   # Check for compiler warnings/remarks
   # grep -v 'feupdateenv ...' ignores a known FDS MPI compiler warning (http://software.intel.com/en-us/forums/showthread.php?t=62806)
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage6f | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 6f warnings:" >> $FIREBOT_DIR/output/warnings
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage6f | grep -v 'feupdateenv is not implemented' | grep -v 'lcilkrts linked' >> $FIREBOT_DIR/output/warnings
      echo "" >> $FIREBOT_DIR/output/warnings
   fi
}

#  ===============================================
#  = Stage 6g - Make SMV pictures (release mode) =
#  ===============================================

make_smv_pictures()
{
   # Run Make SMV Pictures script (release mode)
   cd $SVNROOT/Verification/scripts
   ./Make_SMV_Pictures.sh -r &> $FIREBOT_DIR/output/stage6g
}

check_smv_pictures()
{
   # Scan and report any errors in make SMV pictures process
   cd $FIREBOT_DIR
   if [[ `grep -B 50 -A 50 "Segmentation" -I $FIREBOT_DIR/output/stage6g` == "" ]]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 6g: Make SMV Pictures (Release Mode)"
      grep -B 50 -A 50 "Segmentation" -I $FIREBOT_DIR/output/stage6g > $FIREBOT_DIR/output/stage6g_errors
      ERROR_LOG=$FIREBOT_DIR/output/stage6g_errors
      save_build_status
      email_error_message
   fi
}

#  ================================
#  = Stage 6h - Make FDS pictures =
#  ================================

make_fds_pictures()
{
   # Run Make FDS Pictures script
   cd $SVNROOT/Verification
   ./Make_FDS_Pictures.sh &> $FIREBOT_DIR/output/stage6h
}

check_fds_pictures()
{
   # Scan and report any errors in make FDS pictures process
   cd $FIREBOT_DIR
   if [[ `grep -B 50 -A 50 "Segmentation" -I $FIREBOT_DIR/output/stage6h` == "" ]]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 6h: Make FDS Pictures"
      grep -B 50 -A 50 "Segmentation" -I $FIREBOT_DIR/output/stage6h > $FIREBOT_DIR/output/stage6h_errors
      ERROR_LOG=$FIREBOT_DIR/output/stage6h_errors
      save_build_status
      email_error_message
   fi
}

#  ============================================
#  = Stage 7 - Matlab plotting and statistics =
#  ============================================

run_matlab_plotting()
{
   # Run Matlab plotting script
   cd $SVNROOT/Utilities/Matlab/scripts

   # Replace LaTeX with TeX for Interpreter in plot_style.m
   # This allows displayless automatic Matlab plotting
   # Otherwise Matlab crashes due to a known bug
   sed -i 's/LaTeX/TeX/g' plot_style.m 

   cd $SVNROOT/Utilities/Matlab
   matlab -r "try, disp('Running Matlab Verification script'), FDS_verification_script, catch, disp('Matlab error'), err = lasterror, err.message, err.stack, end, exit" &> $FIREBOT_DIR/output/stage7_verification
   matlab -r "try, disp('Running Matlab Validation script'), FDS_validation_script, catch, disp('Matlab error'), err = lasterror, err.message, err.stack, end, exit" &> $FIREBOT_DIR/output/stage7_validation
}

check_matlab_plotting()
{
   # Scan and report any errors in Matlab scripts
   cd $FIREBOT_DIR
   if [[ `grep -A 50 -E "Matlab error|License checkout failed" $FIREBOT_DIR/output/stage7*` == "" ]]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 7: Matlab plotting and statistics"
      grep -A 50 -E "Matlab error|License checkout failed" $FIREBOT_DIR/output/stage7* > $FIREBOT_DIR/output/stage7_errors
      ERROR_LOG=$FIREBOT_DIR/output/stage7_errors
      save_build_status
      email_error_message
   fi
}

#  ==================================
#  = Stage 8 - Build FDS-SMV Guides =
#  ==================================

make_fds_user_guide()
{
   # Build FDS User Guide
   cd $SVNROOT/Manuals/FDS_User_Guide
   pdflatex -interaction nonstopmode FDS_User_Guide &> $FIREBOT_DIR/output/stage8_fds_user_guide
   bibtex FDS_User_Guide >> $FIREBOT_DIR/output/stage8_fds_user_guide 2>&1
   pdflatex -interaction nonstopmode FDS_User_Guide >> $FIREBOT_DIR/output/stage8_fds_user_guide 2>&1
   pdflatex -interaction nonstopmode FDS_User_Guide >> $FIREBOT_DIR/output/stage8_fds_user_guide 2>&1
}

make_fds_technical_guide()
{
   # Build FDS Technical Guide
   cd $SVNROOT/Manuals/FDS_Technical_Reference_Guide
   pdflatex -interaction nonstopmode FDS_Technical_Reference_Guide &> $FIREBOT_DIR/output/stage8_fds_technical_guide
   bibtex FDS_Technical_Reference_Guide >> $FIREBOT_DIR/output/stage8_fds_technical_guide 2>&1
   pdflatex -interaction nonstopmode FDS_Technical_Reference_Guide >> $FIREBOT_DIR/output/stage8_fds_technical_guide 2>&1
   pdflatex -interaction nonstopmode FDS_Technical_Reference_Guide >> $FIREBOT_DIR/output/stage8_fds_technical_guide 2>&1
}

make_fds_verification_guide()
{
   # Build FDS Verification Guide
   cd $SVNROOT/Manuals/FDS_Verification_Guide
   pdflatex -interaction nonstopmode FDS_Verification_Guide &> $FIREBOT_DIR/output/stage8_fds_verification_guide
   bibtex FDS_Verification_Guide >> $FIREBOT_DIR/output/stage8_fds_verification_guide 2>&1
   pdflatex -interaction nonstopmode FDS_Verification_Guide >> $FIREBOT_DIR/output/stage8_fds_verification_guide 2>&1
   pdflatex -interaction nonstopmode FDS_Verification_Guide >> $FIREBOT_DIR/output/stage8_fds_verification_guide 2>&1
}

make_fds_validation_guide()
{
   # Build FDS Validation Guide
   cd $SVNROOT/Manuals/FDS_Validation_Guide
   pdflatex -interaction nonstopmode FDS_Validation_Guide &> $FIREBOT_DIR/output/stage8_fds_validation_guide
   bibtex FDS_Validation_Guide >> $FIREBOT_DIR/output/stage8_fds_validation_guide 2>&1
   pdflatex -interaction nonstopmode FDS_Validation_Guide >> $FIREBOT_DIR/output/stage8_fds_validation_guide 2>&1
   pdflatex -interaction nonstopmode FDS_Validation_Guide >> $FIREBOT_DIR/output/stage8_fds_validation_guide 2>&1
}

make_fds_configuration_management_plan()
{
   # Build FDS Configuration Management Plan
   cd $SVNROOT/Manuals/FDS_Configuration_Management_Plan
   pdflatex -interaction nonstopmode FDS_Configuration_Management_Plan &> $FIREBOT_DIR/output/stage8_fds_configuration_management_plan
   bibtex FDS_Configuration_Management_Plan >> $FIREBOT_DIR/output/stage8_fds_configuration_management_plan 2>&1
   pdflatex -interaction nonstopmode FDS_Configuration_Management_Plan >> $FIREBOT_DIR/output/stage8_fds_configuration_management_plan 2>&1
   pdflatex -interaction nonstopmode FDS_Configuration_Management_Plan >> $FIREBOT_DIR/output/stage8_fds_configuration_management_plan 2>&1
}

make_smv_user_guide()
{
   # Build SMV User Guide
   cd $SVNROOT/Manuals/SMV_User_Guide
   pdflatex -interaction nonstopmode SMV_User_Guide &> $FIREBOT_DIR/output/stage8_smv_user_guide
   bibtex SMV_User_Guide >> $FIREBOT_DIR/output/stage8_smv_user_guide 2>&1
   pdflatex -interaction nonstopmode SMV_User_Guide >> $FIREBOT_DIR/output/stage8_smv_user_guide 2>&1
   pdflatex -interaction nonstopmode SMV_User_Guide >> $FIREBOT_DIR/output/stage8_smv_user_guide 2>&1
}

make_smv_verification_guide()
{
   # Build SMV Verification Guide
   cd $SVNROOT/Manuals/SMV_Verification_Guide
   pdflatex -interaction nonstopmode SMV_Verification_Guide &> $FIREBOT_DIR/output/stage8_smv_verification_guide
   bibtex SMV_Verification_Guide >> $FIREBOT_DIR/output/stage8_smv_verification_guide 2>&1
   pdflatex -interaction nonstopmode SMV_Verification_Guide >> $FIREBOT_DIR/output/stage8_smv_verification_guide 2>&1
   pdflatex -interaction nonstopmode SMV_Verification_Guide >> $FIREBOT_DIR/output/stage8_smv_verification_guide 2>&1
}

check_all_guides()
{
   # Scan and report any errors in FDS Verification Guide build process
   cd $FIREBOT_DIR
   if [[ `grep "! LaTeX Error:" -I $FIREBOT_DIR/output/stage8*` == "" ]]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 8: FDS-SMV Guides"
      grep "! LaTeX Error:" -I $FIREBOT_DIR/output/stage8* > $FIREBOT_DIR/output/stage8_errors
      ERROR_LOG=$FIREBOT_DIR/output/stage8_errors
      save_build_status
      email_error_message
   fi

   # Check for LaTeX warnings (undefined references or duplicate labels)
   if [[ `grep "undefined|multiply defined|multiply-defined" -I $FIREBOT_DIR/output/stage8*` == "" ]]
   then
      # Continue along
      :
   else
      echo "Stage 8 warnings:" >> $FIREBOT_DIR/output/warnings
      grep "undefined|multiply defined|multiply-defined" -I $FIREBOT_DIR/output/stage8* >> $FIREBOT_DIR/output/warnings
      echo "" >> $FIREBOT_DIR/output/warnings
   fi
}

copy_all_guides_to_website()
{
   # Copy all guides to Blaze status website
   cd $SVNROOT/Manuals
   cp FDS_User_Guide/FDS_User_Guide.pdf \
   FDS_Technical_Reference_Guide/FDS_Technical_Reference_Guide.pdf \
   FDS_Verification_Guide/FDS_Verification_Guide.pdf \
   FDS_Validation_Guide/FDS_Validation_Guide.pdf \
   FDS_Configuration_Management_Plan/FDS_Configuration_Management_Plan.pdf \
   SMV_User_Guide/SMV_User_Guide.pdf \
   SMV_Verification_Guide/SMV_Verification_Guide.pdf \
   /var/www/html/firebot/manuals/
}

#  ==================================================
#  = Build status report - email and save functions =
#  ==================================================

email_success_message()
{
   cd $FIREBOT_DIR
   # Check for compiler warnings
   if [ -e "output/warnings" ]
   then
      # Send email with success message, include compiler warnings
      mail -s "[Firebot] Build success, with compiler warnings. Revision ${SVN_REVISION} passed all build tests." $mailTo < ${FIREBOT_DIR}/output/warnings > /dev/null
   else
      # Send empty email with success message
      mail -s "[Firebot] Build success! Revision ${SVN_REVISION} passed all build tests." $mailTo < /dev/null > /dev/null
   fi
}

email_error_message()
{
   cd $FIREBOT_DIR
   # Check for compiler warnings
   if [ -e "output/warnings" ]
   then
      cat output/warnings >> $ERROR_LOG

      # Send email with failure message and warnings, body of email contains appropriate log file
      mail -s "[Firebot] Build failure, with compiler warnings! Revision ${SVN_REVISION} build failure at ${BUILD_STAGE_FAILURE}." $mailTo < ${ERROR_LOG} > /dev/null
   else
      # Send email with failure message, body of email contains appropriate log file
      mail -s "[Firebot] Build failure! Revision ${SVN_REVISION} build failure at ${BUILD_STAGE_FAILURE}." $mailTo < ${ERROR_LOG} > /dev/null
   fi
   exit
}

save_build_status()
{
   cd $FIREBOT_DIR
   # Save status outcome of build to a text file
   if [[ $BUILD_STAGE_FAILURE != "" ]]
   then
      echo "Revision ${SVN_REVISION} build failure at ${BUILD_STAGE_FAILURE}." > "$FIREBOT_DIR/history/${SVN_REVISION}.txt"
      cat $ERROR_LOG > "$FIREBOT_DIR/history/${SVN_REVISION}_errors.txt"
   else
      if [ -e "output/warnings" ]
         then 
         echo "Revision ${SVN_REVISION} has compiler warnings." > "$FIREBOT_DIR/history/${SVN_REVISION}.txt"
         cat $FIREBOT_DIR/output/warnings > "$FIREBOT_DIR/history/${SVN_REVISION}_warnings.txt"
      else
         echo "Build success! Revision ${SVN_REVISION} passed all build tests." > "$FIREBOT_DIR/history/${SVN_REVISION}.txt"
      fi
   fi
}

#  ============================
#  = Primary script execution =
#  ============================

### Stage 1 ###
clean_svn_repo
do_svn_checkout
check_svn_checkout

### Stage 2a ###
compile_fds_db
check_compile_fds_db

### Stage 2b ###
compile_fds_mpi_db
check_compile_fds_mpi_db

### Stage 3 ###
run_verification_cases_short
check_verification_cases_short

### Stage 4a ###
compile_fds
check_compile_fds

### Stage 4b ###
compile_fds_mpi
check_compile_fds_mpi

### Stage 5 ###
run_verification_cases_long
check_verification_cases_long

### Stage 6a ###
compile_smv_utilities
check_smv_utilities

### Stage 6b ###
compile_smv_test_db
check_compile_smv_test_db

### Stage 6c ###
make_smv_pictures_db
check_smv_pictures_db

### Stage 6d ###
compile_smv_test
check_compile_smv_test

### Stage 6e ###
make_smv_pictures_test
check_smv_pictures_test

### Stage 6f ###
compile_smv
check_compile_smv

### Stage 6g ###
make_smv_pictures
check_smv_pictures

### Stage 6h ###
make_fds_pictures
check_fds_pictures

### Stage 7 ###
run_matlab_plotting
check_matlab_plotting

### Stage 8 ###
make_fds_user_guide
make_fds_technical_guide
make_fds_verification_guide
make_fds_validation_guide
make_fds_configuration_management_plan
make_smv_user_guide
make_smv_verification_guide
check_all_guides
copy_all_guides_to_website

### Success! ###
email_success_message
save_build_status