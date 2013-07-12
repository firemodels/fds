#!/bin/bash

# FDTsbot
# FDTs automatic verification and validation test bot
# Kristopher Overholt
# 8/9/2012

# This script is related to the Fire Dynamics Tools (NUREG 1805)
# validation process. This is independent of the FDS Firebot script.
#
# The FDTsbot script performs the following actions:
#
# 1) Compiles the Utilities/Data_Processing/FDTs.f90 program
# 2) Runs the FDTs program for all FDTs validation cases
# 3) Runs the Utilities/Matlab/FDTs_validation_script.m script
# 4) Compiles the FDTs Validation Guide

SVNROOT=`pwd`/../..

#  =================================================
#  = Stages 1 and 2 - Compile and run FDTs program =
#  =================================================

run_fdts_program()
{
   echo ""
   echo "Running Stages 1 and 2 - Compile and run FDTs program"
   echo ""
   cd $SVNROOT/Validation
   ./Run_FDTs_Cases.sh
}

#  ==========================================
#  = Stage 3 - Run Matlab validation script =
#  ==========================================

# Functions to check for an available Matlab license

run_matlab_license_test()
{
   cd $SVNROOT/Utilities/Firebot
   # Run simple test to see if Matlab license is available
   matlab -r "try, disp('Running Matlab License Check'), catch, disp('License Error'), err = lasterror, err.message, err.stack, end, exit" &> FDTsbot.log
}

scan_matlab_license_test()
{
   cd $SVNROOT/Utilities/Firebot
   # Check for failed license
   if [[ `grep "License checkout failed" FDTsbot.log` == "" ]]
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

   # Replace LaTeX with TeX for Interpreter in plot_style.m
   # This allows displayless automatic Matlab plotting
   # Otherwise Matlab crashes due to a known bug
   # sed -i 's/LaTeX/TeX/g' plot_style.m 

   cd $SVNROOT/Utilities/Matlab
   matlab -r "try, disp('Running Matlab Validation script'), FDTs_validation_script, catch, disp('Error'), err = lasterror, err.message, err.stack, end, exit"

   # Restore LaTeX as plot_style interpreter
   # cd $SVNROOT/Utilities/Matlab/scripts
   # sed -i 's/TeX/LaTeX/g' plot_style.m
}

#  =========================================
#  = Stage 4 - Build FDTs Validation Guide =
#  =========================================

make_fdts_validation_guide()
{
   echo ""
   echo "Running Stage 4 - Build FDTs Validation Guide"
   echo ""
   cd $SVNROOT/Manuals/FDTs_Validation_Guide
   ./make_guide.sh
   cp $SVNROOT/Manuals/FDTs_Validation_Guide/FDTs_Validation_Guide.pdf /var/www/html/firebot/fdts_manuals/
}

#  =================
#  = Final cleanup =
#  =================

cleanup_logs()
{
   cd $SVNROOT/Utilities/Firebot
   rm FDTsbot.log
}

#  ============================
#  = Primary script execution =
#  ============================

### Stages 1 and 2 ###
run_fdts_program

### Stage 3 ###
check_matlab_license_server
run_matlab_plotting

### Stage 4 ###
make_fdts_validation_guide

### Final cleanup ###
cleanup_logs

