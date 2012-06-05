#!/bin/bash -f

# This script runs the FDTs cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/..
export FDTs=$SVNROOT/Utilities/Data_Processing/FDTs
export VDIR=$SVNROOT/Validation

cd $VDIR/CAROLFIRE/FDTs
$FDTs CAROLFIRE_THIEF_Inputs.txt

cd $VDIR/Fleury_Heat_Flux/FDTs
$FDTs Fleury_Heat_Flux_FDTs_Input.txt

cd $VDIR/FM_SNL/FDTs
$FDTs FM_SNL_FDTs_Input.txt

cd $VDIR/LLNL_Enclosure/FDTs
$FDTs LLNL_FDTs_Input.txt

cd $VDIR/NIST_NRC/FDTs
$FDTs NIST_NRC_FDTs_Input.txt

cd $VDIR/Steckler_Compartment/FDTs
$FDTs Steckler_MQH_Inputs.txt

cd $VDIR/UL_NFPRF/FDTs
$FDTs UL_NFPRF_FDTs_Input.txt

cd $VDIR/Vettori_Flat_Ceiling/FDTs
$FDTs Vettori_FDTs_Input.txt

cd $VDIR/VTT/FDTs
$FDTs VTT_FDTs_Input.txt

cd $VDIR/USN_Hangars/FDTs
$FDTs USN_Hangars_FDTs_Input.txt

echo FDTs cases complete