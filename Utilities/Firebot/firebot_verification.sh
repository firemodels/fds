#!/bin/bash

# Firebot
# FDS automatIc veRification and validation tEst bot
# Kristopher Overholt
# 6/22/2012

#  ===================
#  = Input variables =
#  ===================

mailTo="kristopher.overholt@nist.gov"
SVNROOT="/home/kjo/firebot_repo"
FIREBOT_DIR="/home/kjo/firebot"
SVN_REVISION=$1

#  ===========
#  = Warning =
#  ===========

# Warn if running as user other than firebot
if [[ `whoami` == "kjo" ]];
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

#  ====================
#  = Sample functions =
#  ====================
 
# sample_are_there_changes()
# {
#    svn up > $FIREBOT_DIR/tmpsvnup

#    grep "^U" $FIREBOT_DIR/tmpsvnup &> /dev/null
#    if [[ $? == 0 ]]
#    then
#       return 1
#    fi

#    grep "^A" $FIREBOT_DIR/tmpsvnup &> /dev/null
#    if [[ $? == 0 ]]
#    then
#      return 1
#    fi

#    grep "^D" $FIREBOT_DIR/tmpsvnup &> /dev/null
#    if [[ $? == 0 ]]
#    then
#       return 1
#    fi

#    grep "^G" $FIREBOT_DIR/tmpsvnup &> /dev/null
#    if [[ $? == 0 ]]
#    then
#       return 1
#    fi
 
#    return 0
# }

#  ===========================
#  = Stage 1 - SVN functions =
#  ===========================

clean_svn_repo()
{
   # Initialize and start with fresh repo
   # Clean Firebot metafiles
   cd $FIREBOT_DIR
   rm output_*

   # Clean up temporary unversioned and modified versioned repository files
   cd $SVNROOT
   svn revert -Rq *
   svn status --no-ignore | grep '^\?' | sed 's/^\?      //'  | xargs -Ixx rm -rf xx
}

do_svn_checkout()
{
   # If an SVN revision number is specified, then get that revision
   if [[ $SVN_REVISION != "" ]]; then
      echo "Checking out revision r${SVN_REVISION}." > $FIREBOT_DIR/output_stage1
      svn update -r $SVN_REVISION >> $FIREBOT_DIR/output_stage1 2>&1
   # If no SVN revision number is specified, then get the latest revision
   else
      echo "Checking out latest revision." > $FIREBOT_DIR/output_stage1
      svn update >> $FIREBOT_DIR/output_stage1 2>&1
      SVN_REVISION=`tail -n 1 $FIREBOT_DIR/output_stage1 | sed "s/[^0-9]//g"`
   fi
}

check_svn_checkout()
{
   cd $SVNROOT
   # XX grep for errors in $FIREBOT_DIR/output_stage1
}

#  =============================
#  = Stage 2a - Compile FDS DB =
#  =============================

compile_fds_db()
{
   # Clean and compile FDS DB
   cd $SVNROOT/FDS_Compilation/intel_linux_64_db
   make --makefile ../makefile clean &> /dev/null
   ./make_fds.sh &> $FIREBOT_DIR/output_stage2a
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
      ERROR_LOG=$FIREBOT_DIR/output_stage2a
      save_build_status
      email_error_message
   fi

   # Check for compiler warnings
   if [[ `grep warning ${FIREBOT_DIR}/output_stage2a` == "" ]]
   then
      # Continue along
      :
   else
      grep warning ${FIREBOT_DIR}/output_stage2a >> $FIREBOT_DIR/output_compiler_warnings
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
   ./make_fds.sh &> $FIREBOT_DIR/output_stage2b
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
      ERROR_LOG=$FIREBOT_DIR/output_stage2b
      save_build_status
      email_error_message
   fi

   # Check for compiler warnings
   # grep -v 'feupdateenv ...' ignores a known FDS MPI compiler warning (http://software.intel.com/en-us/forums/showthread.php?t=62806)
   if [[ `grep warning ${FIREBOT_DIR}/output_stage2b | grep -v 'feupdateenv is not implemented'` == "" ]]
   then
      # Continue along
      :
   else
      grep warning ${FIREBOT_DIR}/output_stage2b | grep -v 'feupdateenv is not implemented' >> $FIREBOT_DIR/output_compiler_warnings
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
      echo "Waiting for ${JOBS_REMAINING} verification cases to start." >> $FIREBOT_DIR/output_stage3
      sleep 30
   done
}

wait_verification_cases_short_end()
{
   # Scans qstat and waits for verification cases to end
   while [[ `qstat | grep $(whoami)` != '' ]]; do
      JOBS_REMAINING=`qstat | grep $(whoami) | wc -l`
      echo "Waiting for ${JOBS_REMAINING} verification cases to complete." >> $FIREBOT_DIR/output_stage3
      sleep 30
   done
}

run_verification_cases_short()
{
   # Set variables for launching FDS cases on cluster
   cd $SVNROOT/Verification
   export SVNROOT=$SVNROOT
   export FDS=$SVNROOT/FDS_Compilation/intel_linux_64_db/fds_intel_linux_64_db
   export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64_db/fds_mpi_intel_linux_64_db
   export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
   export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh
   export BASEDIR=$SVNROOT/Verification

   #  ========================
   #  = Run all serial cases =
   #  ========================

   # Wait for serial verification cases to start
   ./FDS_Cases.sh &> $FIREBOT_DIR/output_stage3
   wait_verification_cases_short_start

   # Wait some additional time for cases to start
   sleep 30

   # Stop all cases
   export STOPFDS=1
   ./FDS_Cases.sh >> $FIREBOT_DIR/output_stage3 2>&1
   unset STOPFDS

   # Wait for serial verification cases to end
   wait_verification_cases_short_end

   #  =====================
   #  = Run all MPI cases =
   #  =====================

   # Wait for MPI verification cases to start
   ./FDS_MPI_Cases.sh >> $FIREBOT_DIR/output_stage3 2>&1
   wait_verification_cases_short_start

   # Wait some additional time for cases to start
   sleep 30

   # Stop all cases
   export STOPFDS=1
   ./FDS_MPI_Cases.sh >> $FIREBOT_DIR/output_stage3 2>&1
   unset STOPFDS

   # Wait for MPI verification cases to start
   wait_verification_cases_short_end

   # Remove all .stop files from Verification directories (recursively)
   find . -name '*.stop' -exec rm -f {} \;
}

check_verification_cases_short()
{
   # Scan and report any errors in FDS verification cases
   cd $SVNROOT/Verification

   if [[ `grep 'Run aborted' -rI ${FIREBOT_DIR}/output_stage3` == "" ]] && \
      [[ `grep ERROR: -rI *` == "" ]] && \
      [[ `grep 'STOP: Numerical' -rI *` == "" ]] && \
      [[ `grep forrtl -rI *` == "" ]]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 3: FDS Verification Cases"
      
      grep 'Run aborted' -rI $FIREBOT_DIR/output_stage3 >> $FIREBOT_DIR/output_stage3_errors
      grep ERROR: -rI * >> $FIREBOT_DIR/output_stage3_errors
      grep 'STOP: Numerical' -rI * >> $FIREBOT_DIR/output_stage3_errors
      grep forrtl -rI * >> $FIREBOT_DIR/output_stage3_errors
      
      ERROR_LOG=$FIREBOT_DIR/output_stage3_errors
      save_build_status
      email_error_message
   fi
}

#  ==========================
#  = Stage 4a - Compile FDS =
#  ==========================

compile_fds()
{
   # Clean and compile FDS
   cd $SVNROOT/FDS_Compilation/intel_linux_64
   make --makefile ../makefile clean &> /dev/null
   ./make_fds.sh &> $FIREBOT_DIR/output_stage4a
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
      BUILD_STAGE_FAILURE="Stage 4a: FDS Compilation"
      ERROR_LOG=$FIREBOT_DIR/output_stage4a
      save_build_status
      email_error_message
   fi

   # Check for compiler warnings
   if [[ `grep warning ${FIREBOT_DIR}/output_stage4a` == "" ]]
   then
      # Continue along
      :
   else
      grep warning ${FIREBOT_DIR}/output_stage4a >> $FIREBOT_DIR/output_compiler_warnings
   fi
}

#  ==============================
#  = Stage 4a - Compile FDS MPI =
#  ==============================

compile_fds_mpi()
{
   # Clean and compile FDS MPI
   cd $SVNROOT/FDS_Compilation/mpi_intel_linux_64
   make --makefile ../makefile clean &> /dev/null
   ./make_fds.sh &> $FIREBOT_DIR/output_stage4b
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
      BUILD_STAGE_FAILURE="Stage 4b: FDS MPI Compilation"
      ERROR_LOG=$FIREBOT_DIR/output_stage4b
      save_build_status
      email_error_message
   fi

   # Check for compiler warnings
   # grep -v 'feupdateenv ...' ignores a known FDS MPI compiler warning (http://software.intel.com/en-us/forums/showthread.php?t=62806)
   if [[ `grep warning ${FIREBOT_DIR}/output_stage4b | grep -v 'feupdateenv is not implemented'` == "" ]]
   then
      # Continue along
      :
   else
      grep warning ${FIREBOT_DIR}/output_stage4b | grep -v 'feupdateenv is not implemented' >> $FIREBOT_DIR/output_compiler_warnings
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
      echo "Waiting for ${JOBS_REMAINING} verification cases to complete." >> $FIREBOT_DIR/output_stage5
      sleep 60
   done
}

run_verification_cases_long()
{
   # Set variables for launching FDS cases on cluster
   cd $SVNROOT/Verification
   export SVNROOT=$SVNROOT
   export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
   export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64
   export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
   export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh
   export BASEDIR=$SVNROOT/Verification

   # Start running all cases
   ./FDS_Cases.sh &> $FIREBOT_DIR/output_stage5
   ./FDS_MPI_Cases.sh >> $FIREBOT_DIR/output_stage5 2>&1

   # Wait for all verification cases to end
   wait_verification_cases_long_end
}

check_verification_cases_long()
{
   # Scan and report any errors in FDS verification cases
   cd $SVNROOT/Verification

   if [[ `grep 'Run aborted' -rI ${FIREBOT_DIR}/output_stage5` == "" ]] && \
      [[ `grep ERROR: -rI *` == "" ]] && \
      [[ `grep 'STOP: Numerical' -rI *` == "" ]] && \
      [[ `grep forrtl -rI *` == "" ]]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 5: FDS Verification Cases"
      
      grep 'Run aborted' -rI $FIREBOT_DIR/output_stage5 >> $FIREBOT_DIR/output_stage5_errors
      grep ERROR: -rI * >> $FIREBOT_DIR/output_stage5_errors
      grep 'STOP: Numerical' -rI * >> $FIREBOT_DIR/output_stage5_errors
      grep forrtl -rI * >> $FIREBOT_DIR/output_stage5_errors
      
      ERROR_LOG=$FIREBOT_DIR/output_stage5_errors
      save_build_status
      email_error_message
   fi
}

#  =========================================
#  = Stage 6 - Generate Smokeview pictures =
#  =========================================

make_fds_pictures()
{
   # Run Make FDS Pictures script
   cd $SVNROOT/Verification
   ./Make_FDS_Pictures.sh &> $FIREBOT_DIR/output_stage6
}

#  ============================================
#  = Stage 7 - Matlab plotting and statistics =
#  ============================================

run_matlab_plotting()
{
   # Run Matlab plotting script
   cd $SVNROOT/Utilities/Matlab/scripts

   # Replace LaTeX with TeX for Interpreter in plot_style.m
   # This allows displayless autmaotic Matlab plotting
   # Otherwise Matlab crashes due to a known bug
   sed -i 's/LaTeX/TeX/g' plot_style.m 

   cd $SVNROOT/Utilities/Matlab
   matlab -nodisplay -r "FDS_verification_script;quit" &> $FIREBOT_DIR/output_stage7
}

check_matlab_plotting()
{
   cd $SVNROOT/Manuals/FDS_Verification_Guide
   # XX grep matlab output for errors
}

#  ======================================
#  = Stage 8 - Build Verification Guide =
#  ======================================

make_verification_guide()
{
   # Build FDS Verification guide
   cd $SVNROOT/Manuals/FDS_Verification_Guide
   pdflatex -interaction nonstopmode FDS_Verification_Guide &> $FIREBOT_DIR/output_stage8
   bibtex FDS_Verification_Guide >> $FIREBOT_DIR/output_stage8 2>&1
   pdflatex -interaction nonstopmode FDS_Verification_Guide >> $FIREBOT_DIR/output_stage8 2>&1
   pdflatex -interaction nonstopmode FDS_Verification_Guide >> $FIREBOT_DIR/output_stage8 2>&1
}

check_verification_guide()
{
   # Scan and report any errors in FDS Verification Guide build process
   cd $SVNROOT/Manuals/FDS_Verification_Guide
   if [[ `grep "! LaTeX Error:" -I $FIREBOT_DIR/output_stage8` == "" ]]
   then
      # Continue along
      :
   else
      BUILD_STAGE_FAILURE="Stage 8: FDS Verification Guide"
      grep "! LaTeX Error:" -I $FIREBOT_DIR/output_stage8 > $FIREBOT_DIR/output_stage8_errors
      ERROR_LOG=$FIREBOT_DIR/output_stage8_errors
      save_build_status
      email_error_message
   fi
}

#  ==========================
#  = Build status functions =
#  ==========================

email_success_message()
{
   cd $FIREBOT_DIR
   # Check for compiler warnings
   if [ -e "output_compiler_warnings" ]
   then
      # Send email with success message, include compiler warnings
      mail -s "[Firebot] Build success, with compiler warnings. Revision ${SVN_REVISION} passed all build tests." $mailTo < ${FIREBOT_DIR}/output_compiler_warnings > /dev/null
   else
      # Send empty email with success message
      mail -s "[Firebot] Build success! Revision ${SVN_REVISION} passed all build tests." $mailTo < /dev/null > /dev/null
   fi
}

email_error_message()
{
   cd $FIREBOT_DIR
   # Check for compiler warnings
   if [ -e "output_compiler_warnings" ]
   then
      cat output_compiler_warnings >> $ERROR_LOG

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
      if [ -e "output_compiler_warnings" ]
         then 
         echo "Revision ${SVN_REVISION} has compiler warnings." > "$FIREBOT_DIR/history/${SVN_REVISION}.txt"
         cat $FIREBOT_DIR/output_compiler_warnings > "$FIREBOT_DIR/history/${SVN_REVISION}_errors.txt"
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

### Stage 6 ###
# make_fds_pictures

### Stage 7 ###
# run_matlab_plotting
# check_matlab_plotting

### Stage 8 ###
# make_verification_guide
# check_verification_guide

### Success! ###
email_success_message
save_build_status
