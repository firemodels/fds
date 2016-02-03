#!/bin/csh -f

# To include validation data sets to be run automatically by Validationbot,
# be sure to include the data set in ~/FDS-SMV/Validation/Process_All_Output.sh

cd Beyler_Hood; ./Run_All.sh -y; cd ..
cd BRE_Spray; ./Run_All.sh -y; cd ..
cd Bryant_Doorway; ./Run_All.sh -y; cd ..
cd CAROLFIRE; ./Run_All.sh -y; cd ..
cd CHRISTIFIRE; ./Run_All.sh -y; cd ..
cd Cup_Burner; ./Run_All.sh -y; cd ..
cd FAA_Polymers; ./Run_All.sh -y; cd ..
cd Fleury_Heat_Flux; ./Run_All.sh -y; cd ..
cd FM_Parallel_Panels; ./Run_All.sh -y; cd ..
cd Hamins_CH4; ./Run_All.sh -y; cd ..
cd Heskestad_Flame_Height; ./Run_All.sh -y; cd ..
cd LEMTA_Spray; ./Run_All.sh -y; cd ..
cd LLNL_Enclosure; ./Run_All.sh -y; cd ..
cd Moody_Chart; ./Run_All.sh -y; cd ..
cd NIST_He_2009; ./Run_All.sh -y; cd ..
cd NIST_Smoke_Alarms; ./Run_All.sh -y; cd ..
cd NIST_NRC; ./Run_All.sh -y; cd ..
cd NRCC_Facade; ./Run_All.sh -y; cd ..
cd NRL_HAI; ./Run_All.sh -y; cd ..
cd Pool_Fires; ./Run_All.sh -y; cd ..
cd PRISME; ./Run_All.sh -y; cd ..
cd Restivo_Experiment; ./Run_All.sh -y; cd ..
cd Sippola_Aerosol_Deposition; ./Run_All.sh -y; cd ..
cd SP_AST; ./Run_All.sh -y; cd ..
cd Steckler_Compartment; ./Run_All.sh -y; cd ..
cd UL_NIST_Vents; ./Run_All.sh -y; cd ..
cd Ulster_SBI; ./Run_All.sh -y; cd ..
cd UMD_Polymers; ./Run_All.sh -y; cd ..
cd USCG_HAI; ./Run_All.sh -y; cd ..
cd Vettori_Flat_Ceiling; ./Run_All.sh -y; cd ..
cd Vettori_Sloped_Ceiling; ./Run_All.sh -y; cd ..
cd VTT_Sprays; ./Run_All.sh -y; cd ..
cd WTC; ./Run_All.sh -y; cd ..
