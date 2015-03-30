#!/bin/bash

# Firebot variables
FDS_SVNROOT=/home/jgs/FDS-SMV
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

function usage {
echo "build_2waycode.sh [ -q queue_name -r revision_number -s -u svn_username -v max_validation_processes -y ]"
echo "Runs testing script for 2waycode"
echo ""
echo "Options"
echo "-r - revision_number - run cases using a specific SVN revision number"
echo "     default: (none, latest SVN HEAD)"
echo ""
exit
}

# Update repository
SVN_REVISION=''
while getopts 'hq:r:su:v:y' OPTION
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
esac
done
shift $(($OPTIND-1))

if [[ $SVN_REVISION = "" ]]; then
   cd $FDS_SVNROOT/FDS_Source
   svn update >> $OUTPUT_DIR/stage1 2>&1
   cd $TWOWAY_DIR/source
   svn update >> $OUTPUT_DIR/stage1_twowaycode 2>&1
else
   cd $FDS_SVNROOT/FDS_Source
   svn update -r $SVN_REVISION >> $OUTPUT_DIR/stage1 2>&1
   cd $TWOWAY_DIR/source
   svn update -r $SVN_REVISION >> $OUTPUT_DIR/stage1_twowaycode 2>&1
   echo "At revision ${SVN_REVISION}." 
fi

SVN_REVISION=`tail -n 1 $OUTPUT_DIR/stage1 | sed "s/[^0-9]//g"`
echo $SVN_REVISION

# Print the FDS revision number on User Guide
#cd $TWOWAY_DIR
#sed -i "s:.*SVN Repository Revision.*:SVN Repository Revision ${SVN_REVISION}:" fds2ftmi_user_guide.tex

# Print the FDS revision number on python scripts
#cd $FIREBOT_DIR
#sed -i "s:.*SVN=.*:SVN='${SVN_REVISION}':" generate_plots.py

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
rm twowaycode_user_guide.tex
svn up -r $SVN_REVISION twowaycode_user_guide.tex

# Revert the FDS revision number on python scripts
cd $FIREBOT_DIR
rm generate_plots.py
svn up -r $SVN_REVISION generate_plots.py

exit