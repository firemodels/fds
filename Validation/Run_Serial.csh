#!/bin/csh -f

# To include validation data sets to be run automatically by Validationbot,
# be sure to include the data set in ~/FDS-SMV/Validation/Process_All_Output.sh

cd Beyler_Hood; mkdir Current_Results; ./Run_All.sh; cd ..
cd BRE_Spray; mkdir Current_Results; ./Run_All.sh; cd ..
cd Bryant_Doorway; mkdir Current_Results; ./Run_All.sh; cd ..
cd CAROLFIRE; mkdir Current_Results; ./Run_All.sh; cd ..
cd CHRISTIFIRE; mkdir Current_Results; ./Run_All.sh; cd ..
cd Cup_Burner; mkdir Current_Results; ./Run_All.sh; cd ..
cd FAA_Polymers; mkdir Current_Results; ./Run_All.sh; cd ..
cd Fleury_Heat_Flux; mkdir Current_Results; ./Run_All.sh; cd ..
cd FM_Parallel_Panels; mkdir Current_Results; ./Run_All.sh; cd ..
cd Hamins_CH4; mkdir Current_Results; ./Run_All.sh; cd ..
cd Heskestad_Flame_Height; mkdir Current_Results; ./Run_All.sh; cd ..
cd LEMTA_Spray; mkdir Current_Results; ./Run_All.sh; cd ..
cd LLNL_Enclosure; mkdir Current_Results; ./Run_All.sh; cd ..
cd Moody_Chart; mkdir Current_Results; ./Run_All.sh; cd ..
cd NIST_He_2009; mkdir Current_Results; ./Run_All.sh; cd ..
cd NIST_Smoke_Alarms; mkdir Current_Results; ./Run_All.sh; cd ..
cd NIST_NRC; mkdir Current_Results; ./Run_All.sh; cd ..
cd NRCC_Facade; mkdir Current_Results; ./Run_All.sh; cd ..
cd NRL_HAI; mkdir Current_Results; ./Run_All.sh; cd ..
cd Pool_Fires; mkdir Current_Results; ./Run_All.sh; cd ..
cd PRISME; mkdir Current_Results; ./Run_All.sh; cd ..
cd Restivo_Experiment; mkdir Current_Results; ./Run_All.sh; cd ..
cd Sippola_Aerosol_Deposition; mkdir Current_Results; ./Run_All.sh; cd ..
cd SP_AST; mkdir Current_Results; ./Run_All.sh; cd ..
cd Steckler_Compartment; mkdir Current_Results; ./Run_All.sh; cd ..
cd UL_Commodity; mkdir Current_Results; ./Run_All.sh; cd ..
cd UL_NIST_Vents; mkdir Current_Results; ./Run_All.sh; cd ..
cd Ulster_SBI; mkdir Current_Results; ./Run_All.sh; cd ..
cd USCG_HAI; mkdir Current_Results; ./Run_All.sh; cd ..
cd Vettori_Flat_Ceiling; mkdir Current_Results; ./Run_All.sh; cd ..
cd Vettori_Sloped_Ceiling; mkdir Current_Results; ./Run_All.sh; cd ..
cd VTT_Sprays; mkdir Current_Results; ./Run_All.sh; cd ..
cd WTC; mkdir Current_Results; ./Run_All.sh; cd ..
