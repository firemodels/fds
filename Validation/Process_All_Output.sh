#!/bin/bash

# This script contains a list of validation sets whose results are to be copied to the 'out' repository.
# The script checks to see that every case within each validation set has run successfully.
# The script also checks to see if older versions of the results might accidentally overwrite newer versions.

# This function scans the '_git.txt' files within a given directory and converts the git hash to a Unix time stamp.

GET_UNIX_TIME()
{
  gitrevisions=$TEMPDIR/gitrevisions.$$
  cat $1/*git.txt 2> /dev/null | sort -u > $gitrevisions
  gitrev=`head -1 $gitrevisions`
  gitdate2=0
  if [ "$gitrev" != "" ] ; then
    gitrevshort=`echo $gitrev | awk -F - '{print $(NF-1)}' | sed 's/^g//'`
    gitdate2=`git show -s --format=%at $gitrevshort | head -1 | awk '{print $1}'`
  fi
  rm $gitrevisions
}

# This function takes a validation set and determines whether or not it should be processed.

PROCESS()
{
  case=$1
  casedir=$FIREMODELS/out/$case
  GET_UNIX_TIME $casedir
  out_date=$gitdate2
  curdir=`pwd`
  casedir=$case/Current_Results
  GET_UNIX_TIME $casedir
  new_date=$gitdate2
  cd $case
  nout=`ls -l Current_Results/*.out |& grep -v cannot | wc -l`
  nfds=`ls -l Current_Results/*.fds |& grep -v cannot | wc -l`
  ncfds=`ls -l Current_Results/*cat.fds |& grep -v cannot | wc -l`
  if [ $ncfds -gt 0 ] ; then
    nfds=$ncfds
  fi
  nsuccess=`tail Current_Results/*.out |& grep successfully | wc -l`
  status="***error: $case cases not run"
  if [ $nfds -gt 0 ] && [ $nfds -gt $nout ]; then
    status="***error: some $case cases did not run or are not complete"
  elif [ $out_date \> $new_date ]; then
    status="***error: existing output is newer than the cases being processed"
  else
    if [ $nout -gt 0 ] && [ $nout -gt $nsuccess ]; then
      status="some $case cases failed"
    else
      if [ $nout -gt 0 ] ; then
      status="processing output"
      ./Process_Output.sh
      fi
    fi
  fi
  if [ $nfds -gt 0 ]; then
    echo "$case: cases=$nfds finished=$nout successful=$nsuccess status=$status"
  else
    echo "$case: No cases run"
  fi
  cd $curdir
}

# Define some directories.

TEMPDIR=$HOME/temp
if [ ! -d $TEMPDIR ]; then
   mkdir $TEMPDIR
fi

FIREMODELS=../..

# There should exist a line entry for every directory in the Validation directory of the fds repository.
# If the case is under development, simply comment out the line below.

PROCESS Aalto_Woods
PROCESS Arup_Tunnel
PROCESS ATF_Corridors
PROCESS Atmospheric_Dispersion
PROCESS Backward_Facing_Step
PROCESS Beyler_Hood
PROCESS BGC_GRI_LNG_Fires
PROCESS Bittern_Sprinkler_Experiments
PROCESS Bouchair_Solar_Chimney
PROCESS BRE_Spray
PROCESS Bryant_Doorway
PROCESS CAROLFIRE
PROCESS Casara_Arts_Ribbed_Channel
PROCESS Convection
PROCESS Crown_Fires
PROCESS CSIRO_Grassland_Fires
PROCESS CSTB_Tunnel
PROCESS Cup_Burner
PROCESS DelCo_Trainers
PROCESS DoJ_HAI_Pool_Fires
PROCESS Droplet_Evaporation
PROCESS Edinburgh_Vegetation_Drag
PROCESS FAA_Cargo_Compartments
PROCESS FAA_Polymers
PROCESS Fleury_Heat_Flux
PROCESS FM_Burner
PROCESS FM_FPRF_Datacenter
PROCESS FM_Parallel_Panels
PROCESS FM_SNL
PROCESS FM_Vertical_Wall_Flames
PROCESS Frankman_Vegetation
PROCESS Hamins_Gas_Burners
PROCESS Harrison_Spill_Plumes
PROCESS Heated_Channel_Flow
PROCESS Heskestad_Flame_Height
PROCESS Insulation_Materials
PROCESS JH_FRA
PROCESS Juelich_SETCOM
PROCESS LEMTA_Spray
PROCESS LEMTA_UGent_Pool_Fires
PROCESS LLNL_Enclosure
PROCESS LNG_Dispersion
PROCESS Loughborough_Jet_Fires
PROCESS McCaffrey_Plume
PROCESS Memorial_Tunnel
PROCESS Montoir_LNG_Fires
PROCESS Moody_Chart
PROCESS MPI_Scaling_Tests
PROCESS NBS_Multi-Room
PROCESS NIST_Composite_Beam
PROCESS NIST_Deposition_Gauge
PROCESS NIST_Douglas_Firs
PROCESS NIST_E119_Compartment
PROCESS NIST_FSE_2008
PROCESS NIST_He_2009
PROCESS NIST_NRC
PROCESS NIST_NRC_Corner_Effects
PROCESS NIST_NRC_OLIVE-Fire
PROCESS NIST_NRC_Parallel_Panels
PROCESS NIST_NRC_Transient_Combustibles
PROCESS NIST_Polymers
PROCESS NIST_Pool_Fires
PROCESS NIST_RSE_1994
PROCESS NIST_RSE_2007
PROCESS NIST_Smoke_Alarms
PROCESS NIST_Structure_Separation
PROCESS NIST_Vent_Study
PROCESS NRCC_Facade
PROCESS NRCC_Smoke_Tower
PROCESS NRL_HAI
PROCESS OMP_Scaling_Tests
PROCESS Phoenix_LNG_Fires
PROCESS Pool_Fires
PROCESS PRISME
PROCESS Purdue_Flames
PROCESS Ranz_Marshall
PROCESS Restivo_Experiment
PROCESS Sandia_Methane_Burner
PROCESS Sandia_Plumes
PROCESS Scaling_Pyrolysis
PROCESS Shell_LNG_Fireballs
PROCESS Sippola_Aerosol_Deposition
PROCESS Smyth_Slot_Burner
PROCESS SP_AST
PROCESS SP_Wood_Cribs
PROCESS Steckler_Compartment
PROCESS SWJTU_Tunnels
PROCESS Theobald_Hose_Stream
PROCESS Turbulent_Jet
PROCESS TUS_Facade
PROCESS UL_NFPRF
PROCESS UL_NIJ_Houses
PROCESS UL_NIST_Vents
PROCESS Ulster_SBI
PROCESS UMD_Burning_Rate_Emulator
PROCESS UMD_Line_Burner
PROCESS UMD_Polymers
PROCESS UMD_SBI
PROCESS USCG_HAI
PROCESS USFS_Catchpole
PROCESS USFS_Corsica
PROCESS USN_Hangars
PROCESS UWO_Wind_Tunnel
PROCESS Vettori_Flat_Ceiling
PROCESS Vettori_Sloped_Ceiling
PROCESS VTT
PROCESS VTT_Sprays
PROCESS Waterloo_Methanol
PROCESS WTC
PROCESS Wu_Bakar_Tunnels
