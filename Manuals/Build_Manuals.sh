#!/bin/bash

# This script is handy for locally building all the manuals at once for a release bundle

export manuals_dir=`pwd`

cd $manuals_dir/FDS_Config_Management_Plan
echo "Building FDS Configuration Management Plan at: " $manuals_dir/FDS_Config_Management_Plan
./make_guide.sh
cp FDS_Config_Management_Plan.pdf $manuals_dir/.

cd $manuals_dir/FDS_User_Guide
echo "Building FDS User Guide at: " $manuals_dir/FDS_User_Guide
./make_guide.sh
cp FDS_User_Guide.pdf $manuals_dir/.

cd $manuals_dir/FDS_Technical_Reference_Guide
echo "Building FDS Technical Reference Guide at: " $manuals_dir/FDS_Technical_Reference_Guide
./make_guide.sh
cp FDS_Technical_Reference_Guide.pdf $manuals_dir/.

cd $manuals_dir/FDS_Verification_Guide
echo "Building FDS Verification Guide at: " $manuals_dir/FDS_Verification_Guide
./make_guide.sh
cp FDS_Verification_Guide.pdf $manuals_dir/.

cd $manuals_dir/FDS_Validation_Guide
echo "Building FDS Validation Guide at: " $manuals_dir/FDS_Validation_Guide
./make_guide.sh
cp FDS_Validation_Guide.pdf $manuals_dir/.

