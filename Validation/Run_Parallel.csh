#!/bin/csh -f

# To include validation data sets to be run automatically by Validationbot,
# be sure to include the data set in ~/FDS-SMV/Validation/Process_All_Output.sh

cd Arup_Tunnel;  ./Run_All.sh -y; cd ..
cd ATF_Corridors;  ./Run_All.sh -y; cd ..
cd Backward_Facing_Step;  ./Run_All.sh -y; cd ..
cd Bouchair_Solar_Chimney;  ./Run_All.sh -y; cd ..
cd DelCo_Trainers;  ./Run_All.sh -y; cd ..
cd FAA_Cargo_Compartments;  ./Run_All.sh -y; cd ..
cd FM_FPRF_Datacenter;  ./Run_All.sh -y; cd ..
cd FM_SNL;  ./Run_All.sh -y; cd ..
cd FM_Vertical_Wall_Flames;  ./Run_All.sh -y; cd ..
cd Harrison_Spill_Plumes;  ./Run_All.sh -y; cd ..
cd Heated_Channel_Flow;  ./Run_All.sh -y; cd ..
cd LNG_Dispersion;  ./Run_All.sh -y; cd ..
cd McCaffrey_Plume;  ./Run_All.sh -y; cd ..
cd NBS_Multi-Room;  ./Run_All.sh -y; cd ..
cd NIST_Douglas_Firs;  ./Run_All.sh -y; cd ..
cd NIST_FSE_2008;  ./Run_All.sh -y; cd ..
cd NIST_NRC_Corner_Effects;  ./Run_All.sh -y; cd ..
cd NIST_RSE_1994;  ./Run_All.sh -y; cd ..
cd NIST_RSE_2007;  ./Run_All.sh -y; cd ..
cd NIST_Vent_Study;  ./Run_All.sh -y; cd ..
cd Purdue_Flames;  ./Run_All.sh -y; cd ..
cd Sandia_Plumes;  ./Run_All.sh -y; cd ..
cd Smyth_Slot_Burner;  ./Run_All.sh -y; cd ..
cd Turbulent_Jet;  ./Run_All.sh -y; cd ..
cd UL_NFPRF;  ./Run_All.sh -y; cd ..
cd UL_NIJ_Houses;  ./Run_All.sh -y; cd ..
cd USN_Hangars;  ./Run_All.sh -y; cd ..
cd VTT;  ./Run_All.sh -y; cd ..
cd Waterloo_Methanol;  ./Run_All.sh -y; cd ..
