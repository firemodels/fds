#!/bin/bash -f

# This script runs all of the FDTs cases

export SVNROOT=`pwd`/..
export FDTs=$SVNROOT/Utilities/Data_Processing/FDTs
export VERIFICATION_DIR=$SVNROOT/Verification
export VALIDATION_DIR=$SVNROOT/Validation

# First, compile latest version of FDTs.f90

cd $SVNROOT/Utilities/Data_Processing
ifort FDTs.f90 -o FDTs

# Run all FDTs verification cases

cd $VERIFICATION_DIR/FDTs
$FDTs Verification.txt

# Run all FDTs validation cases

cd $VALIDATION_DIR/ATF_Corridors/FDTs
$FDTs ATF_Corridors_FDTs_Inputs.txt

cd $VALIDATION_DIR/CAROLFIRE/FDTs
$FDTs CAROLFIRE_FDTs_Inputs.txt

cd $VALIDATION_DIR/Fleury_Heat_Flux/FDTs
$FDTs Fleury_Heat_Flux_FDTs_Inputs.txt

cd $VALIDATION_DIR/FM_SNL/FDTs
$FDTs FM_SNL_FDTs_Inputs.txt

cd $VALIDATION_DIR/LLNL_Enclosure/FDTs
$FDTs LLNL_FDTs_Inputs.txt

cd $VALIDATION_DIR/NBS_Multi-Room/FDTs
$FDTs NBS_Multi-Room_FDTs_Inputs.txt

cd $VALIDATION_DIR/NIST_NRC/FDTs
$FDTs NIST_NRC_FDTs_Inputs.txt

cd $VALIDATION_DIR/NIST_Dunes_2000/FDTs
$FDTs NIST_Dunes_2000_FDTs_Inputs.txt

cd $VALIDATION_DIR/SP_AST/FDTs
$FDTs SP_AST_FDTs_Inputs.txt

cd $VALIDATION_DIR/Steckler_Compartment/FDTs
$FDTs Steckler_FDTs_Inputs.txt

cd $VALIDATION_DIR/UL_NFPRF/FDTs
$FDTs UL_NFPRF_FDTs_Inputs.txt

cd $VALIDATION_DIR/USN_Hangars/FDTs
$FDTs USN_Hangars_FDTs_Inputs.txt

cd $VALIDATION_DIR/Vettori_Flat_Ceiling/FDTs
$FDTs Vettori_FDTs_Inputs.txt

cd $VALIDATION_DIR/VTT/FDTs
$FDTs VTT_FDTs_Inputs.txt

cd $VALIDATION_DIR/WTC/FDTs
$FDTs WTC_FDTs_Inputs.txt

echo FDTs cases complete
