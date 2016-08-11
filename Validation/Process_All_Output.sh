#!/bin/bash -f

if [ "$FDSSMV" == "" ] ; then
  export FDSSMV=~/FDS-SMVgitclean
fi
export VDIR=$FDSSMV/Validation

# This list of active validation data sets is used by Validationbot
# to automatically run validation cases on a regular basis.

# There should exist a line entry for every directory under Validation.
# If the case is under development, simply comment out the line.

$VDIR/Arup_Tunnel/FDS_Output_Files/Process_Output.csh
$VDIR/ATF_Corridors/FDS_Output_Files/Process_Output.csh
$VDIR/Backward_Facing_Step/FDS_Output_Files/Process_Output.csh
$VDIR/Beyler_Hood/FDS_Output_Files/Process_Output.csh
$VDIR/BRE_Spray/FDS_Output_Files/Process_Output.csh
$VDIR/Bryant_Doorway/FDS_Output_Files/Process_Output.csh
$VDIR/CAROLFIRE/FDS_Output_Files/Process_Output.csh
$VDIR/CHRISTIFIRE/FDS_Output_Files/Process_Output.csh
$VDIR/CSIRO_Grassland_Fires/FDS_Output_Files/Process_Output.csh
$VDIR/Cup_Burner/FDS_Output_Files/Process_Output.csh
#$VDIR/DelCo_Trainers/FDS_Output_Files/Process_Output.csh
$VDIR/FAA_Cargo_Compartments/FDS_Output_Files/Process_Output.csh
$VDIR/FAA_Polymers/FDS_Output_Files/Process_Output.csh
$VDIR/Fleury_Heat_Flux/FDS_Output_Files/Process_Output.csh
$VDIR/FM_Parallel_Panels/FDS_Output_Files/Process_Output.csh
$VDIR/FM_SNL/FDS_Output_Files/Process_Output.csh
$VDIR/Hamins_Gas_Burners/FDS_Output_Files/Process_Output.csh
$VDIR/Harrison_Spill_Plumes/FDS_Output_Files/Process_Output.csh
$VDIR/Heskestad_Flame_Height/FDS_Output_Files/Process_Output.csh
$VDIR/LEMTA_Spray/FDS_Output_Files/Process_Output.csh
$VDIR/LLNL_Enclosure/FDS_Output_Files/Process_Output.csh
$VDIR/McCaffrey_Plume/FDS_Output_Files/Process_Output.csh
$VDIR/Moody_Chart/FDS_Output_Files/Process_Output.csh
$VDIR/MPI_Scaling_Tests/FDS_Output_Files/Process_Output.csh
$VDIR/NBS_Multi-Room/FDS_Output_Files/Process_Output.csh
$VDIR/NIST_Douglas_Firs/FDS_Output_Files/Process_Output.csh
$VDIR/NIST_FSE_2008/FDS_Output_Files/Process_Output.csh
$VDIR/NIST_He_2009/FDS_Output_Files/Process_Output.csh
$VDIR/NIST_NRC/FDS_Output_Files/Process_Output.csh
$VDIR/NIST_RSE_1994/FDS_Output_Files/Process_Output.csh
$VDIR/NIST_Smoke_Alarms/FDS_Output_Files/Process_Output.csh
$VDIR/NRCC_Facade/FDS_Output_Files/Process_Output.csh
$VDIR/NRCC_Smoke_Tower/FDS_Output_Files/Process_Output.csh
$VDIR/NRL_HAI/FDS_Output_Files/Process_Output.csh
$VDIR/Pool_Fires/FDS_Output_Files/Process_Output.csh
$VDIR/PRISME/FDS_Output_Files/Process_Output.csh
$VDIR/Purdue_Flames/FDS_Output_Files/Process_Output.csh
$VDIR/Restivo_Experiment/FDS_Output_Files/Process_Output.csh
$VDIR/Sandia_Plumes/FDS_Output_Files/Process_Output.csh
$VDIR/Sippola_Aerosol_Deposition/FDS_Output_Files/Process_Output.csh
$VDIR/Smyth_Slot_Burner/FDS_Output_Files/Process_Output.csh
$VDIR/SP_AST/FDS_Output_Files/Process_Output.csh
$VDIR/Steckler_Compartment/FDS_Output_Files/Process_Output.csh
$VDIR/Turbulent_Jet/FDS_Output_Files/Process_Output.csh
$VDIR/UL_NFPRF/FDS_Output_Files/Process_Output.csh
$VDIR/UL_NIST_Vents/FDS_Output_Files/Process_Output.csh
$VDIR/Ulster_SBI/FDS_Output_Files/Process_Output.csh
$VDIR/UMD_Polymers/FDS_Output_Files/Process_Output.csh
$VDIR/UMD_Line_Burner/FDS_Output_Files/Process_Output.csh
$VDIR/USCG_HAI/FDS_Output_Files/Process_Output.csh
$VDIR/USN_Hangars/FDS_Output_Files/Process_Output.csh
$VDIR/Vettori_Flat_Ceiling/FDS_Output_Files/Process_Output.csh
$VDIR/Vettori_Sloped_Ceiling/FDS_Output_Files/Process_Output.csh
$VDIR/VTT/FDS_Output_Files/Process_Output.csh
$VDIR/VTT_Sprays/FDS_Output_Files/Process_Output.csh
$VDIR/WTC/FDS_Output_Files/Process_Output.csh
