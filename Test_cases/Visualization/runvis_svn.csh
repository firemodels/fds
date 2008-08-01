#!/bin/csh -f

# This script, runvis_svn.csh, is a template under svn control.
# Before using, this template should first be copied and customized 
# for each set test cases to be run.

# specify location of the fds5 executables
setenv FDS ~/bin/fds5_linux
setenv FDSMPI ~/bin/fds5_mpi_linux

# specify location of repository root
setenv SVNROOT ~/FDS-SMV

# Option:
# Uncomment the setenv line below to stop all FDS jobs running via this script.
# setenv STOPFDS

#  1.  To run this script, first copy runall_svn.csh to runall.csh  
#      (only when runall_svn.csh changes)
#  2.  define the FDS and/or FDSMPI environment variables to point
#      to the version of fds you want to run.
#  3.  Change hostnames in each RUNFDS (or RUNFDS_MPI) command to point to free 
#      cluster nodes
#  4.  Run the script in the directory CONTAINING directories listed 
#      on the various RUNFDS/RUNFDSMPI command lines.

# VVVVVVVVVVV Do not change these line VVVVVVVVVVVVVV
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
set RUNFDSMPI=$SVNROOT/Scripts/runfdsmpi.csh
setenv BASEDIR `pwd`
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# syntax of RUNFDS
# $RUNFDS  directory case host

# To use with MPI define the LAMNODES environment variable to point
# the the symbolic names of the cluster nodes to be used.
#
# mpi example
# -----------
# setenv LAMNODES n0 n0 n0 n0
# $RUNFDSMPI Demonstrations/2Room_ranch ranch_00 fire72 &

# VVVVVVVVVVVV Change lines below to point to cases to be run VVVVVVVVVVVVV

# plume cases
$RUNFDS . plume5a fire64 &
$RUNFDS . plume5b fire64 &
$RUNFDS . plume5c fire65 &
$RUNFDS . pplume5 fire65 &
setenv LAMNODES n0 n0 n0 n0
$RUNFDSMPI . thouse5 fire73 &

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
