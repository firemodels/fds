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
   sed -i 's/LaTeX/TeX/g' plot_style.m 

   cd $SVNROOT/Utilities/Matlab
   matlab -r "try, disp('Running Matlab Validation script'), FDTs_validation_script, catch, disp('Error'), err = lasterror, err.message, err.stack, end, exit"

   # Restore LaTeX as plot_style interpreter
   cd $SVNROOT/Utilities/Matlab/scripts
   sed -i 's/TeX/LaTeX/g' plot_style.m
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

#  ============================
#  = Primary script execution =
#  ============================

### Stages 1 and 2 ###
run_fdts_program

### Stage 3 ###
run_matlab_plotting

### Stage 4 ###
make_fdts_validation_guide

