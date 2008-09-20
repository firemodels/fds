#!/bin/csh -f
#
# This script performs a "ls -l" for each FDS .out file
# in order to quickly determine which cases have finished
# 

# specify location of repository root
set SVNROOT=~/FDS-SMV

set getstatus=$SVNROOT/Utilities/Scripts/display_status.csh

cd $SVNROOT/Training
date
# demonstration cases
$getstatus Demonstrations/2Room_Ranch ranch_00
$getstatus Demonstrations/2Room_Ranch ranch_01
$getstatus Demonstrations/2Room_Ranch ranch_02
$getstatus Demonstrations/2Room_Ranch ranch_03
$getstatus Demonstrations/2Room_Ranch ranch_04
# MCFRS cases
$getstatus MCFRS/MCFRS_Flashover MCFRS_Flashover_00
#$getstatus MCFRS/MCFRS_Flashover MCFRS_Flashover_01
#$getstatus MCFRS/MCFRS_Flashover MCFRS_Flashover_02
#$getstatus MCFRS/MCFRS_Flashover MCFRS_Flashover_03
$getstatus MCFRS/MCFRS_Ranch MCFRS_Ranch_00
# MFRI  cases
$getstatus MFRI/MFRI_Training_Tower MFRI_Training_Tower_00
$getstatus MFRI/MFRI_Training_Tower MFRI_Training_Tower_01
$getstatus MFRI/MFRI_Training_Tower MFRI_Training_Tower_02
$getstatus MFRI/MFRI_Training_Tower MFRI_Training_Tower_03
$getstatus MFRI/MFRI_Training_Tower MFRI_Training_Tower_04
$getstatus MFRI/MFRI_Training_Tower MFRI_Training_Tower_05
$getstatus MFRI/MFRI_Training_Tower MFRI_Training_Tower_06
$getstatus MFRI/MFRI_Training_Tower MFRI_Training_Tower_07
