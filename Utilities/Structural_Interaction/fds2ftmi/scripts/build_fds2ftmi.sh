#!/bin/bash

# Firebot variables
FDS_SVNROOT=/home/jgs/FDS-SMV
FIREBOT_DIR=$FDS_SVNROOT/Utilities/Structural_Interaction/fds2ftmi/scripts
FDS2FTMI_DIR=$FDS_SVNROOT/Utilities/Structural_Interaction/fds2ftmi
OUTPUT_DIR=/home/jgs/FDS-SMV/Utilities/Structural_Interaction/fds2ftmi/scripts/output
ERROR_LOG=$OUTPUT_DIR/errors
WARNING_LOG=$OUTPUT_DIR/warnings

# Clean outputs
cd $FIREBOT_DIR
rm -rf output/

# Create output dir
mkdir output

function usage {
echo "firebot.sh [ -q queue_name -r revision_number -s -u svn_username -v max_validation_processes -y ]"
echo "Runs Firebot V&V testing script"
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
else
   cd $FDS_SVNROOT/FDS_Source
   svn update -r $SVN_REVISION >> $OUTPUT_DIR/stage1 2>&1
   echo "At revision ${SVN_REVISION}." >> $OUTPUT_DIR/stage1 2>&1
fi

SVN_REVISION=`tail -n 1 $OUTPUT_DIR/stage1 | sed "s/[^0-9]//g"`
echo $SVN_REVISION

# Print the FDS revision number on User Guide
cd $FDS2FTMI_DIR
sed -i "s:.*SVN Repository Revision.*:SVN Repository Revision ${SVN_REVISION}:" fds2ftmi_user_guide.tex

# Print the FDS revision number on python scripts
cd $FIREBOT_DIR
sed -i "s:.*SVN=.*:SVN='${SVN_REVISION}':" generate_plots.py

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

# Compile debug version of fds2ftmi
compile_fds2ftmi_db()
{
   # Clean and compile FDS2ftmi debug
   cd $FDS_SVNROOT/Utilities/Structural_Interaction/fds2ftmi/intel_linux_64_db
   make -f ../makefile clean &> /dev/null
   ./make_fds2ftmi.sh &> $OUTPUT_DIR/stage2a_ftmi
}

check_compile_fds2ftmi_db()
{
   # Check for errors in FDS debug compilation
   cd $FDS_SVNROOT/Utilities/Structural_Interaction/fds2ftmi/intel_linux_64_db
   if [ -e "fds2ftmi_linux_64_db" ]
   then
      stage2a_ftmi_success=true
   else
      echo "Errors from Stage 2a_ftmi - Compile and inspect FDS2ftmi debug:" >> $ERROR_LOG
      cat ${FIREBOT_DIR}/output/stage2a_ftmi >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2a_ftmi` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 2a_ftmi - Compile and inspect FDS2ftmi debug:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage2a_ftmi >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

# Compile release version of fds2ftmi
compile_fds2ftmi()
{
   # Clean and compile FDS2ftmi debug
   cd $FDS_SVNROOT/Utilities/Structural_Interaction/fds2ftmi/intel_linux_64
   make -f ../makefile clean &> /dev/null
   ./make_fds2ftmi.sh &> $OUTPUT_DIR/stage4a_ftmi
}

check_compile_fds2ftmi()
{
   # Check for errors in FDS debug compilation
   cd $FDS_SVNROOT/Utilities/Structural_Interaction/fds2ftmi/intel_linux_64
   if [ -e "fds2ftmi_linux_64" ]
   then
      stage4a_ftmi_success=true
   else
      echo "Errors from Stage 4a_ftmi - Compile FDS2ftmi release:" >> $ERROR_LOG
      cat ${FIREBOT_DIR}/output/stage4a_ftmi >> $ERROR_LOG
      echo "" >> $ERROR_LOG
   fi

   # Check for compiler warnings/remarks
   if [[ `grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage4a_ftmi | grep -v 'performing multi-file optimizations' | grep -v 'generating object file'` == "" ]]
   then
      # Continue along
      :
   else
      echo "Warnings from Stage 4a_ftmi - Compile FDS2ftmi release:" >> $WARNING_LOG
      grep -A 5 -E 'warning|remark' ${FIREBOT_DIR}/output/stage4a_ftmi | grep -v 'performing multi-file optimizations' | grep -v 'generating object file' >> $WARNING_LOG
      echo "" >> $WARNING_LOG
   fi
}

# Functions to check for an available Ansys license

run_ansys_license_test()
{
   # Run simple test to see if ansys license is available
   cd $FDS_SVNROOT/Utilities/Structural_Interaction/fds2ftmi/scripts
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
      TIME_LIMIT_STAGE="7"
      check_time_limit
      # Wait 5 minutes until retry
      sleep 300
      check_ansys_license_server
   fi
}

check_ansys_license_server()
{
   run_ansys_license_test
   scan_ansys_license_test
}

# Compile and check release and debug versions of FDS
compile_fds_db
check_compile_fds_db
compile_fds
check_compile_fds

# Compile and check debug and release versions
compile_fds2ftmi_db
check_compile_fds2ftmi_db
compile_fds2ftmi
check_compile_fds2ftmi

# Check Ansys license status
run_ansys_license_test
scan_ansys_license_test

# Run Verification Cases
# simple_panel_hot
export OMP_NUM_THREADS=2
cd $FDS2FTMI_DIR/examples/simple_panel_hot
./simple_panel_hot.sh
# h_profile
cd $FDS2FTMI_DIR/examples/h_profile
./h_profile.sh

# Run Verification plots
cd $FDS2FTMI_DIR/scripts
python generate_plots.py

# Build Guide
cd $FDS2FTMI_DIR
export TEXINPUTS=$TEXINPUTS:.:../../../Manuals/LaTeX_Style_Files/
pdflatex fds2ftmi_user_guide.tex
bibtex fds2ftmi_user_guide.tex
pdflatex fds2ftmi_user_guide.tex

exit