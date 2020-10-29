#!/bin/bash

# To include validation data sets to be run automatically by Validationbot,
# be sure to include the data set in ~/FDS-SMV/Validation/Process_All_Output.sh

OPTIONS="$* -y"

CURDIR=`pwd`

cd $CURDIR/Arup_Tunnel;                ./Run_All.sh $OPTIONS
cd $CURDIR/ATF_Corridors;              ./Run_All.sh $OPTIONS
cd $CURDIR/Atmospheric_Dispersion;     ./Run_All.sh $OPTIONS
cd $CURDIR/Backward_Facing_Step;       ./Run_All.sh $OPTIONS
cd $CURDIR/Bouchair_Solar_Chimney;     ./Run_All.sh $OPTIONS
cd $CURDIR/Crown_Fires;                ./Run_All.sh $OPTIONS
cd $CURDIR/CSIRO_Grassland_Fires;      ./Run_All.sh $OPTIONS
cd $CURDIR/CSTB_Tunnel;                ./Run_All.sh $OPTIONS
cd $CURDIR/DelCo_Trainers;             ./Run_All.sh $OPTIONS
cd $CURDIR/FAA_Cargo_Compartments;     ./Run_All.sh $OPTIONS
cd $CURDIR/FM_Burner;                  ./Run_All.sh $OPTIONS
cd $CURDIR/FM_FPRF_Datacenter;         ./Run_All.sh $OPTIONS
cd $CURDIR/FM_SNL;                     ./Run_All.sh $OPTIONS
cd $CURDIR/FM_Vertical_Wall_Flames;    ./Run_All.sh $OPTIONS
cd $CURDIR/Harrison_Spill_Plumes;      ./Run_All.sh $OPTIONS
cd $CURDIR/Heated_Channel_Flow;        ./Run_All.sh $OPTIONS
cd $CURDIR/Juelich_SETCOM;             ./Run_All.sh $OPTIONS
cd $CURDIR/LNG_Dispersion;             ./Run_All.sh $OPTIONS
cd $CURDIR/McCaffrey_Plume;            ./Run_All.sh $OPTIONS
cd $CURDIR/Moody_Chart;                ./Run_All.sh $OPTIONS
cd $CURDIR/NBS_Multi-Room;             ./Run_All.sh $OPTIONS
cd $CURDIR/NIST_Composite_Beam;        ./Run_All.sh $OPTIONS
cd $CURDIR/NIST_Deposition_Gauge;      ./Run_All.sh $OPTIONS
cd $CURDIR/NIST_Douglas_Firs;          ./Run_All.sh $OPTIONS
cd $CURDIR/NIST_E119_Compartment;      ./Run_All.sh $OPTIONS
cd $CURDIR/NIST_FSE_2008;              ./Run_All.sh $OPTIONS
cd $CURDIR/NIST_NRC_Corner_Effects;    ./Run_All.sh $OPTIONS
cd $CURDIR/NIST_NRC_Parallel_Panels;   ./Run_All.sh $OPTIONS
cd $CURDIR/NIST_Pool_Fires;            ./Run_All.sh $OPTIONS
cd $CURDIR/NIST_RSE_1994;              ./Run_All.sh $OPTIONS
cd $CURDIR/NIST_RSE_2007;              ./Run_All.sh $OPTIONS
cd $CURDIR/NIST_Vent_Study;            ./Run_All.sh $OPTIONS
cd $CURDIR/NRCC_Smoke_Tower;           ./Run_All.sh $OPTIONS
cd $CURDIR/Pool_Fires;                 ./Run_All.sh $OPTIONS
cd $CURDIR/Purdue_Flames;              ./Run_All.sh $OPTIONS
cd $CURDIR/Sandia_Plumes;              ./Run_All.sh $OPTIONS
cd $CURDIR/Sippola_Aerosol_Deposition; ./Run_All.sh $OPTIONS
cd $CURDIR/Smyth_Slot_Burner;          ./Run_All.sh $OPTIONS
cd $CURDIR/SWJTU_Tunnels;              ./Run_All.sh $OPTIONS
cd $CURDIR/Turbulent_Jet;              ./Run_All.sh $OPTIONS
cd $CURDIR/UL_NFPRF;                   ./Run_All.sh $OPTIONS
cd $CURDIR/UL_NIJ_Houses;              ./Run_All.sh $OPTIONS
cd $CURDIR/UMD_Line_Burner;            ./Run_All.sh $OPTIONS
cd $CURDIR/USFS_Catchpole;             ./Run_All.sh $OPTIONS
cd $CURDIR/USFS_Corsica;               ./Run_All.sh $OPTIONS
cd $CURDIR/USN_Hangars;                ./Run_All.sh $OPTIONS
cd $CURDIR/UWO_Wind_Tunnel;            ./Run_All.sh $OPTIONS
cd $CURDIR/VTT;                        ./Run_All.sh $OPTIONS
cd $CURDIR/Waterloo_Methanol;          ./Run_All.sh $OPTIONS
