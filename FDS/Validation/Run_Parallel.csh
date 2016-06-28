#!/bin/csh -f

# To include validation data sets to be run automatically by Validationbot,
# be sure to include the data set in ~/FDS-SMV/Validation/Process_All_Output.sh

cd Arup_Tunnel;  ./Run_All.sh -y; cd ..
cd ATF_Corridors;  ./Run_All.sh -y; cd ..
cd Backward_Facing_Step;  ./Run_All.sh -y; cd ..
cd CSIRO_Grassland_Fires;  ./Run_All.sh -y; cd ..
cd FAA_Cargo_Compartments;  ./Run_All.sh -y; cd ..
cd FM_FPRF_Datacenter;  ./Run_All.sh -y; cd ..
cd FM_SNL;  ./Run_All.sh -y; cd ..
cd Harrison_Spill_Plumes;  ./Run_All.sh -y; cd ..
cd McCaffrey_Plume;  ./Run_All.sh -y; cd ..
cd NBS_Multi-Room;  ./Run_All.sh -y; cd ..
cd NIST_Douglas_Firs;  ./Run_All.sh -y; cd ..
cd NIST_FSE_2008;  ./Run_All.sh -y; cd ..
cd NIST_RSE_1994;  ./Run_All.sh -y; cd ..
cd Purdue_Flames;  ./Run_All.sh -y; cd ..
cd Sandia_Plumes;  ./Run_All.sh -y; cd ..
cd Smyth_Slot_Burner;  ./Run_All.sh -y; cd ..
cd Turbulent_Jet;  ./Run_All.sh -y; cd ..
cd UL_NFPRF;  ./Run_All.sh -y; cd ..
cd UMD_Line_Burner;  ./Run_All.sh -y; cd ..
cd USN_Hangars;  ./Run_All.sh -y; cd ..
cd VTT;  ./Run_All.sh -y; cd ..
