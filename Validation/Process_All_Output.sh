#!/bin/bash

PROCESS()
{
  case=$1
  curdir=`pwd`
  cd $case
  nout=`ls -l Current_Results/*.out |& grep -v cannot | wc -l`
  nfds=`ls -l Current_Results/*.fds |& grep -v cannot | wc -l`
  nsuccess=`tail Current_Results/*.out |& grep successfully | wc -l`
  status="***error: $case cases not run"
  if [ $nfds -gt 0 ] && [ $nfds -gt $nout ]; then
    status="***error: some $case cases did not run or are not complete"
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

# This list of active validation data sets is used by Validationbot
# to automatically run validation cases on a regular basis.

# There should exist a line entry for every directory under Validation.
# If the case is under development, simply comment out the line.

PROCESS Arup_Tunnel
PROCESS ATF_Corridors
PROCESS Backward_Facing_Step
PROCESS Beyler_Hood
PROCESS Bouchair_Solar_Chimney
PROCESS BRE_Spray
PROCESS Bryant_Doorway
PROCESS CAROLFIRE
PROCESS CHRISTIFIRE
PROCESS CSIRO_Grassland_Fires
PROCESS Cup_Burner
PROCESS DelCo_Trainers
PROCESS FAA_Cargo_Compartments
PROCESS FAA_Polymers
PROCESS Fleury_Heat_Flux
PROCESS FM_FPRF_Datacenter
PROCESS FM_Parallel_Panels
PROCESS FM_SNL
PROCESS FM_Vertical_Wall_Flames
PROCESS Hamins_Gas_Burners
PROCESS Harrison_Spill_Plumes
PROCESS Heated_Channel_Flow
PROCESS Heskestad_Flame_Height
PROCESS LEMTA_Spray
PROCESS LLNL_Enclosure
PROCESS LNG_Dispersion
PROCESS McCaffrey_Plume
PROCESS Moody_Chart
PROCESS MPI_Scaling_Tests
PROCESS NBS_Multi-Room
PROCESS NIST_Douglas_Firs
PROCESS NIST_FSE_2008
PROCESS NIST_He_2009
PROCESS NIST_NRC
PROCESS NIST_NRC_Corner_Effects
PROCESS NIST_RSE_1994
PROCESS NIST_RSE_2007
PROCESS NIST_Smoke_Alarms
PROCESS NIST_Vent_Study
PROCESS NRCC_Facade
PROCESS NRCC_Smoke_Tower
PROCESS NRL_HAI
PROCESS Pool_Fires
PROCESS PRISME
PROCESS Purdue_Flames
PROCESS Restivo_Experiment
PROCESS Sandia_Plumes
PROCESS Sippola_Aerosol_Deposition
PROCESS Smyth_Slot_Burner
PROCESS SP_AST
PROCESS Steckler_Compartment
PROCESS Turbulent_Jet
PROCESS UL_NFPRF
PROCESS UL_NIJ_Houses
PROCESS UL_NIST_Vents
PROCESS Ulster_SBI
PROCESS UMD_Line_Burner
PROCESS UMD_Polymers
PROCESS USCG_HAI
PROCESS USN_Hangars
PROCESS Vettori_Flat_Ceiling
PROCESS Vettori_Sloped_Ceiling
PROCESS VTT
PROCESS VTT_Sprays
PROCESS WTC
