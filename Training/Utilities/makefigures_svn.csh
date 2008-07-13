#!/bin/csh -f

# This script, makefigures_svn.csh, is a template under svn control.
# It runs smokeview for each test case generating desired figures
# to be included in the Trainer report.  This script
# should be modified as test cases are added or removed.

# specify location of the smokeview executable
setenv SMV5 ~/bin/smv5_intel

#  1.  To use this script, first copy makefigures_svn.csh to makefigures.csh  
#      (only when makefigures_svn.csh changes)
#  2.  define the SMV5 environment variable to point to the version 
#      of smokeview you want to run.
#  3.  Run script from this directory (repository_root/bin)

# only edit hostnames below

set RUNSMV=./runsmv.csh
setenv JOBDIR `pwd`/..
setenv REPORTFIGDIR `pwd`/../../Manuals/FDS_SMV_Training_Guide/FIGURES/SCRIPT

# syntax of RUNFDS
# $RUNFDS  directory case host

# demonstration cases
#$RUNSMV Demonstrations/2Room_Ranch ranch_00 
#$RUNSMV Demonstrations/2Room_Ranch ranch_01
#$RUNSMV Demonstrations/2Room_Ranch ranch_02
#$RUNSMV Demonstrations/2Room_Ranch ranch_03
#$RUNSMV Demonstrations/2Room_Ranch ranch_04
# MCFRS cases
#$RUNSMV MCFRS/MCFRS_Flashover MCFRS_Flashover_00
$RUNSMV MCFRS/MCFRS_Flashover MCFRS_Flashover_01
#$RUNSMV MCFRS/MCFRS_Flashover MCFRS_Flashover_02
#$RUNSMV MCFRS/MCFRS_Flashover MCFRS_Flashover_03
#$RUNSMV MCFRS/MCFRS_Ranch MCFRS_Ranch_00
# MFRI  cases
#$RUNSMV MFRI/MFRI_Training_Tower MFRI_Training_Tower_00
