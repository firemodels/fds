#!/bin/csh -f

# To include validation data sets to be run automatically by Validationbot,
# be sure to include the data set in ~/FDS-SMV/Validation/Process_All_Output.sh

cd Arup_Tunnel; mkdir Current_Results; ./Run_All.sh; cd ..
cd ATF_Corridors; mkdir Current_Results; ./Run_All.sh; cd ..
cd Backward_Facing_Step; mkdir Current_Results; ./Run_All.sh; cd ..
cd FAA_Cargo_Compartments; mkdir Current_Results; ./Run_All.sh; cd ..
cd FM_SNL; mkdir Current_Results; ./Run_All.sh; cd ..
cd Harrison_Spill_Plumes; mkdir Current_Results; ./Run_All.sh; cd ..
cd McCaffrey_Plume; mkdir Current_Results; ./Run_All.sh; cd ..
cd NBS_Multi-Room; mkdir Current_Results; ./Run_All.sh; cd ..
cd NIST_FSE_2008; mkdir Current_Results; ./Run_All.sh; cd ..
cd NIST_RSE_1994; mkdir Current_Results; ./Run_All.sh; cd ..
cd Purdue_Flames; mkdir Current_Results; ./Run_All.sh; cd ..
cd Sandia_Plumes; mkdir Current_Results; ./Run_All.sh; cd ..
cd Smyth_Slot_Burner; mkdir Current_Results; ./Run_All.sh; cd ..
cd Turbulent_Jet; mkdir Current_Results; ./Run_All.sh; cd ..
cd UL_NFPRF; mkdir Current_Results; ./Run_All.sh; cd ..
cd USN_Hangars; mkdir Current_Results; ./Run_All.sh; cd ..
cd VTT; mkdir Current_Results; ./Run_All.sh; cd ..
