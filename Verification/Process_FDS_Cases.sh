#!/bin/bash -f

# This script runs a few Fortran programs that process some of the output
# of the FDS Verification Cases.

export SVNROOT=`pwd`/..
cd NS_Analytical_Solution
$SVNROOT/Utilities/Data_Processing/ns2d
cd ..
cd Radiation
$SVNROOT/Utilities/Data_Processing/radiation_box
$SVNROOT/Utilities/Data_Processing/radiation_plane_layer
$SVNROOT/Utilities/Data_Processing/wall_internal_radiation
cd ..


