#!/bin/bash -f

# This script runs all of the Correlation cases

export SVNROOT=`pwd`/..
export correlations=$SVNROOT/Utilities/Data_Processing/correlations
export VERIFICATION_DIR=$SVNROOT/Verification
export VALIDATION_DIR=$SVNROOT/Validation

# First, compile latest version of correlations.f90

cd $SVNROOT/Utilities/Data_Processing
ifort correlations.f90 -o correlations

# Run all Correlation verification cases

cd $VERIFICATION_DIR/Correlations
$correlations Verification.txt

# Run all Correlation validation cases

cd $VALIDATION_DIR/ATF_Corridors/Correlations
$correlations ATF_Corridors_Correlation_Inputs.txt

cd $VALIDATION_DIR/CAROLFIRE/Correlations
$correlations CAROLFIRE_Correlation_Inputs.txt

cd $VALIDATION_DIR/Fleury_Heat_Flux/Correlations
$correlations Fleury_Heat_Flux_Correlation_Inputs.txt

cd $VALIDATION_DIR/FM_SNL/Correlations
$correlations FM_SNL_Correlation_Inputs.txt

cd $VALIDATION_DIR/iBMB/Correlations
$correlations iBMB_Correlation_Inputs.txt

cd $VALIDATION_DIR/LLNL_Enclosure/Correlations
$correlations LLNL_Correlation_Inputs.txt

cd $VALIDATION_DIR/NBS_Multi-Room/Correlations
$correlations NBS_Multi-Room_Correlation_Inputs.txt

cd $VALIDATION_DIR/NIST_NRC/Correlations
$correlations NIST_NRC_Correlation_Inputs.txt

cd $VALIDATION_DIR/NIST_Smoke_Alarms/Correlations
$correlations NIST_Smoke_Alarms_Correlation_Inputs.txt

cd $VALIDATION_DIR/SP_AST/Correlations
$correlations SP_AST_Correlation_Inputs.txt

cd $VALIDATION_DIR/Steckler_Compartment/Correlations
$correlations Steckler_Correlation_Inputs.txt

cd $VALIDATION_DIR/UL_NFPRF/Correlations
$correlations UL_NFPRF_Correlation_Inputs.txt

cd $VALIDATION_DIR/USN_Hangars/Correlations
$correlations USN_Hangars_Correlation_Inputs.txt

cd $VALIDATION_DIR/Vettori_Flat_Ceiling/Correlations
$correlations Vettori_Correlation_Inputs.txt

cd $VALIDATION_DIR/VTT/Correlations
$correlations VTT_Correlation_Inputs.txt

cd $VALIDATION_DIR/WTC/Correlations
$correlations WTC_Correlation_Inputs.txt

echo Correlation cases complete
