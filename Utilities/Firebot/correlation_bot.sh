#!/bin/bash

# correlation_bot.sh
# Correlations automatic verification and validation test bot
# Kristopher Overholt
# 8/9/2012

# This script is related to the correlation validation process.
# This is independent of the FDS Firebot script.
#
# The correlation_bot script performs the following actions:
#
# 1) Compiles the Utilities/Data_Processing/correlations.f90 program
# 2) Runs the correlations program for all correlation verification and
#    validation cases
# 3) Runs the Utilities/Matlab/Correlation_validation_script.m script
# 4) Compiles the Correlation Guide

SVNROOT=`pwd`/../..

#  =========================================================
#  = Stages 1 and 2 - Compile and run correlations program =
#  =========================================================

run_correlations()
{
   echo ""
   echo "Running Stages 1 and 2 - Compile and run correlations program"
   echo ""
   cd $SVNROOT/Validation
   ./Run_Correlation_Cases.sh
}

#  ==========================================
#  = Stage 3 - Run Matlab validation script =
#  ==========================================

# Functions to check for an available Matlab license

run_matlab_license_test()
{
   cd $SVNROOT/Utilities/Firebot
   # Run simple test to see if Matlab license is available
   matlab -r "try, disp('Running Matlab License Check'), catch, disp('License Error'), err = lasterror, err.message, err.stack, end, exit" &> correlation_bot.log
}

scan_matlab_license_test()
{
   cd $SVNROOT/Utilities/Firebot
   # Check for failed license
   if [[ `grep "License checkout failed" correlation_bot.log` == "" ]]
   then
      # Continue along
      :
   else
      # Wait 5 minutes until retry
      echo "No Matlab licenses available; will retry in 5 minutes."
      sleep 300
      check_matlab_license_server
   fi
}

check_matlab_license_server()
{
   run_matlab_license_test
   scan_matlab_license_test
}

run_matlab_plotting()
{
   echo ""
   echo "Running Stage 3 - Run Matlab validation script"
   echo ""

   # Run Matlab plotting script
   cd $SVNROOT/Utilities/Matlab/scripts

   cd $SVNROOT/Utilities/Matlab
   matlab -r "try, disp('Running Matlab Validation script'), Correlation_validation_script, catch, disp('Error'), err = lasterror, err.message, err.stack, end, exit"
}

#  =====================================
#  = Stage 4 - Build Correlation Guide =
#  =====================================

make_correlation_guide()
{
   echo ""
   echo "Running Stage 4 - Build Correlation Guide"
   echo ""
   cd $SVNROOT/Manuals/Correlation_Guide
   ./make_guide.sh
#  cp $SVNROOT/Manuals/Correlation_Guide/Correlation_Guide.pdf /var/www/html/firebot/correlation_guide/
}

#  =================
#  = Final cleanup =
#  =================

cleanup_logs()
{
   cd $SVNROOT/Utilities/Firebot
   rm correlation_bot.log
}

#  ============================
#  = Primary script execution =
#  ============================

### Stages 1 and 2 ###
run_correlations

### Stage 3 ###
check_matlab_license_server
run_matlab_plotting

### Stage 4 ###
make_correlation_guide

### Final cleanup ###
cleanup_logs

