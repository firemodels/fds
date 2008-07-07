#!/bin/csh -f

# This script, runall_svn.csh, is a template under svn control.
# It contains entries for each test case to be run.  This script
# should be modified as test cases are added or removed.

# specify location of the fds5 executable
setenv FDS5 ~/bin/fds5_intel

# To run a case in parallel, put the following line 
# (filled out correctly) before the case to be run and
#   reset back to the regular FDS5 value afterwards
#
# setenv FDS5 "mpirun n0 n0 n0 n0 ~/bin/fds5_mpi_intel"

# Uncomment the setenv line below to stop all FDS jobs running 
# via this script.
# setenv STOPFDS

#  1.  To run this script, first copy runall_svn.csh to runall.csh  
#      (only when runall_svn.csh changes)
#  2.  define the FDS5 environment variable to point to the version 
#      of fds you want to run.
#  3.  Change hostnames in each RUNFDS command to point to free 
#      cluster nodes
#  4.  Run script from this directory (repository_root/bin)

# only edit hostnames below

set RUNFDS=./runfds.csh
setenv JOBDIR `pwd`/..

# syntax of RUNFDS
# $RUNFDS  directory case host

# demonstration cases
$RUNFDS Demonstrations/2Room_Ranch ranch_00 fire72 &
$RUNFDS Demonstrations/2Room_Ranch ranch_01 fire72 &
$RUNFDS Demonstrations/2Room_Ranch ranch_02 fire72 &
$RUNFDS Demonstrations/2Room_Ranch ranch_03 fire72 &
$RUNFDS Demonstrations/2Room_Ranch ranch_04 fire73 &
# MCFRS cases
$RUNFDS MCFRS/MCFRS_Flashover MCFRS_Flashover_00 fire73 &
$RUNFDS MCFRS/MCFRS_Ranch MCFRS_Ranch_00 fire73 &
# MFRI  cases
$RUNFDS MFRI/MFRI_Training_Tower MFRI_Training_Tower_00 fire74 &
