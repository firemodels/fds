#!/bin/csh -f

# This script, runall_svn.csh, is a template under svn control.
# It contains entries for each test case to be run.  This script
# should be modified as test cases are added or removed.

# specify location of the fds5 executable
setenv FDS5 ~/bin/fds5_intel

#  1.  To run this script, first copy runall_svn.csh to runall.csh  
#      (only when runall_svn.csh changes)
#  2.  define the FDS5 environment variable to point to the version 
#      of fds you want to run.
#  3.  Change hostnames in each RUNFDS command to point to free 
#      cluster nodes
#  4.  Run script from this directory (repository_root/bin)

set RUNFDS=./runfds.csh
setenv JOBDIR `pwd`/..

# syntax of RUNFDS
# $RUNFDS  directory case host

# demontration cases
$RUNFDS demonstrations/2room_ranch ranch_00 fire72 &
# MCFRS cases
$RUNFDS MCFRS/MCFRS_flashover mcfrs_flashover_00 fire72 &
$RUNFDS MCFRS/MCFRS_ranch mcfrs_ranch_00 fire72 &
# MFRI  cases
$RUNFDS MFRI/training_tower mfri fire72 &
$RUNFDS MFRI/training_tower mfri_tower_00a fire72 &
$RUNFDS MFRI/training_tower mfri_tower_00b fire72 &
