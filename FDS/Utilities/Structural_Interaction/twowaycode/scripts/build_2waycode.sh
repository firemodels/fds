#!/bin/bash

# Firebot variables
FIREBOT_DIR=$FDS_SVNROOT/Utilities/Structural_Interaction/twowaycode/scripts
TWOWAY_DIR=$FDS_SVNROOT/Utilities/Structural_Interaction/twowaycode
OUTPUT_DIR=/home/jgs/FDS-SMV/Utilities/Structural_Interaction/twowaycode/scripts/output
ERROR_LOG=$OUTPUT_DIR/errors
WARNING_LOG=$OUTPUT_DIR/warnings

# Clean outputs
cd $FIREBOT_DIR
rm -rf output

# Create output dir
mkdir output

# Clean previous results
cd $TWOWAY_DIR/examples/simply_beam
rm *.csv 

# Get Git Hash
GIT_HASH=$(shell git describe --long --dirty)
echo %GIT_HASH%

# Print the FDS revision number on User Guide
cd $TWOWAY_DIR
sed -i "s:.*Git Hash.*:\\path{%GIT_HASH%}:" twowaycode_user_guide.tex

# Print the FDS revision number on python scripts
cd $FIREBOT_DIR
sed -i "s:.*GIT=.*:GIT='%GIT_HASH%':" generate_plots.py

compile_fds_db()
{
   # Clean and compile FDS debug
   cd $FDS_SVNROOT/FDS_Compilation/intel_linux_64_db
   make -f ../makefile clean &> /dev/null
   ./make_fds.sh &> $OUTPUT_DIR/stage2a
}

check_compile_fds_db()
{
   # Check for errors in FDS debug compilation
   cd $FDS_SVNROOT/FDS_Compilation/intel_linux_64_db
   if [ -e "fds_intel_linux_64_db" ]
   then
      stage2a_success=true
   else
      echo "Errors from Stage 2a - Compile and inspect FDS debug:" >> $ERROR_LOG
      cat $OUTPUT_DIR/stage2a >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2a` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 2a - Compile and inspect FDS debug:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2a >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}
compile_fds()
{
   # Clean and compile FDS
   cd $FDS_SVNROOT/FDS_Compilation/intel_linux_64
   make -f ../makefile clean &> /dev/null
   ./make_fds.sh &> $OUTPUT_DIR/stage4a
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
      cat $OUTPUT_DIR/stage4a >> $ERROR_LOG
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

# Compile debug version of twowaycode
compile_twowaycode_db()
{
   # Clean and compile twowaycode debug
   cd $FDS_SVNROOT/Utilities/Structural_Interaction/twowaycode/intel_linux_64_db
   make -f ../makefile clean &> /dev/null
   ./make_twowaycode.sh &> $OUTPUT_DIR/stage2a_twowaycode
}

check_compile_twowaycode_db()
{
   # Check for errors in FDS debug compilation
   cd $FDS_SVNROOT/Utilities/Structural_Interaction/twowaycode/intel_linux_64_db
   if [ -e "twowaycode_linux_64_db" ]
   then
      stage2a_twowaycode_success=true
   else
      echo "Errors from Stage 2a_twowaycode - Compile and inspect twowaycode debug:" >> $ERROR_LOG
      cat ${FIREBOT_DIR}/output/stage2a_twowaycode >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2a_twowaycode` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 2a_twowaycode - Compile and inspect twowaycode debug:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2a_twowaycode >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

# Compile release version of twowaycode
compile_twowaycode()
{
   # Clean and compile twowaycode debug
   cd $FDS_SVNROOT/Utilities/Structural_Interaction/twowaycode/intel_linux_64
   make -f ../makefile clean &> /dev/null
   ./make_twowaycode.sh &> $OUTPUT_DIR/stage4a_twowaycode
}

check_compile_twowaycode()
{
   # Check for errors in FDS debug compilation
   cd $FDS_SVNROOT/Utilities/Structural_Interaction/twowaycode/intel_linux_64
   if [ -e "twowaycode_linux_64" ]
   then
      stage4a_twowaycode_success=true
   else
      echo "Errors from Stage 4a_twowaycode - Compile twowaycode release:" >> $ERROR_LOG
      cat ${FIREBOT_DIR}/output/stage4a_twowaycode >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage4a_twowaycode | grep -v 'performing multi-file optimizations' | grep -v 'generating object file'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 4a_twowaycode - Compile twowaycode release:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage4a_twowaycode | grep -v 'performing multi-file optimizations' | grep -v 'generating object file' >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

# Functions to check for an available Ansys license

run_ansys_license_test()
{
   # Run simple test to see if ansys license is available
   cd $FDS_SVNROOT/Utilities/Structural_Interaction/twowaycode/scripts
   ansys150 -j lic_test <lic_test.ans> $OUTPUT_DIR/stage7_ansys_license
   rm lic_test.db
   rm lic_test.err
   rm lic_test.log
}

scan_ansys_license_test()
{
   # Check for failed license
   if [[ `grep "ERROR - ANSYS license not available." $OUTPUT_DIR/stage7_ansys_license` == "" ]]
   then
      # Continue along
      :
   else
      # Wait 5 minutes until retry
      sleep 300
      echo "ANSYS license not found!"
      echo "Sleeping for 5 minutes..."
      wait_ansys_lic=$(($wait_ansys_lic + 1 ))    
      check_ansys_license_server
   fi
}

check_ansys_license_server()
{      
   rm $OUTPUT_DIR/stage7_ansys_license
   echo $wait_ansys_lic
   if [ $wait_ansys_lic == 6 ]; then
      echo "ANSYS license error!"
      echo "Exiting..."
      exit
   else
      # Continue along
      :
   fi
   run_ansys_license_test
   scan_ansys_license_test
}

# Compile and check release and debug versions of FDS
compile_fds_db
check_compile_fds_db
compile_fds
check_compile_fds

# Compile and check debug and release versions
compile_twowaycode_db
check_compile_twowaycode_db
compile_twowaycode
check_compile_twowaycode

# Check Ansys license status
wait_ansys_lic=0
run_ansys_license_test
scan_ansys_license_test

# Run Verification Cases
# simply_beam
export OMP_NUM_THREADS=2
cd $TWOWAY_DIR/examples/simply_beam
./simply_beam.sh


# Run Verification plots
cd $TWOWAY_DIR/scripts
python generate_plots.py

# Build Guide
cd $TWOWAY_DIR
export TEXINPUTS=$TEXINPUTS:.:../../../Manuals/LaTeX_Style_Files/
pdflatex twowaycode_user_guide.tex
bibtex twowaycode_user_guide
pdflatex twowaycode_user_guide.tex

# Revert the FDS revision number on User Guide
cd $TWOWAY_DIR
git checkout -- twowaycode_user_guide.tex

# Revert the FDS revision number on python scripts
cd $FIREBOT_DIR
git checkout -- generate_plots.py

exit