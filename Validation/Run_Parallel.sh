#!/bin/bash

# To include validation data sets to be run automatically by Validationbot,
# be sure to include the data set in ~/FDS-SMV/Validation/Process_All_Output.sh

OPTIONS="$* -y"

cd Arup_Tunnel;  ./Run_All.sh $OPTIONS; cd ..
cd ATF_Corridors;  ./Run_All.sh $OPTIONS; cd ..
cd Atmospheric_Dispersion;  ./Run_All.sh $OPTIONS; cd ..
cd Backward_Facing_Step;  ./Run_All.sh $OPTIONS; cd ..
cd Bouchair_Solar_Chimney;  ./Run_All.sh $OPTIONS; cd ..
cd CSIRO_Grassland_Fires;  ./Run_All.sh $OPTIONS; cd ..
cd DelCo_Trainers;  ./Run_All.sh $OPTIONS; cd ..
cd FAA_Cargo_Compartments;  ./Run_All.sh $OPTIONS; cd ..
cd FM_FPRF_Datacenter;  ./Run_All.sh $OPTIONS; cd ..
cd FM_SNL;  ./Run_All.sh $OPTIONS; cd ..
cd FM_Vertical_Wall_Flames;  ./Run_All.sh $OPTIONS; cd ..
cd Harrison_Spill_Plumes;  ./Run_All.sh $OPTIONS; cd ..
cd Heated_Channel_Flow;  ./Run_All.sh $OPTIONS; cd ..
cd LNG_Dispersion;  ./Run_All.sh $OPTIONS; cd ..
cd McCaffrey_Plume;  ./Run_All.sh $OPTIONS; cd ..
cd Moody_Chart; ./Run_All.sh $OPTIONS; cd ..
cd NBS_Multi-Room;  ./Run_All.sh $OPTIONS; cd ..
cd NIST_Composite_Beam;  ./Run_All.sh $OPTIONS; cd ..
cd NIST_Douglas_Firs;  ./Run_All.sh $OPTIONS; cd ..
cd NIST_FSE_2008;  ./Run_All.sh $OPTIONS; cd ..
cd NIST_NRC_Corner_Effects;  ./Run_All.sh $OPTIONS; cd ..
cd NIST_RSE_1994;  ./Run_All.sh $OPTIONS; cd ..
cd NIST_RSE_2007;  ./Run_All.sh $OPTIONS; cd ..
cd NIST_Vent_Study;  ./Run_All.sh $OPTIONS; cd ..
cd NRCC_Smoke_Tower;  ./Run_All.sh $OPTIONS; cd ..
cd Pool_Fires;  ./Run_All.sh $OPTIONS; cd ..
cd Purdue_Flames;  ./Run_All.sh $OPTIONS; cd ..
cd Sandia_Plumes;  ./Run_All.sh $OPTIONS; cd ..
cd Smyth_Slot_Burner;  ./Run_All.sh $OPTIONS; cd ..
cd Turbulent_Jet;  ./Run_All.sh $OPTIONS; cd ..
cd UL_NFPRF;  ./Run_All.sh $OPTIONS; cd ..
cd UL_NIJ_Houses;  ./Run_All.sh $OPTIONS; cd ..
cd UMD_Line_Burner;  ./Run_All.sh $OPTIONS; cd ..
cd USFS_Catchpole;  ./Run_All.sh $OPTIONS; cd ..
cd USN_Hangars;  ./Run_All.sh $OPTIONS; cd ..
cd UWO_Wind_Tunnel;  ./Run_All.sh $OPTIONS; cd ..
cd VTT;  ./Run_All.sh $OPTIONS; cd ..
cd Waterloo_Methanol;  ./Run_All.sh $OPTIONS; cd ..
