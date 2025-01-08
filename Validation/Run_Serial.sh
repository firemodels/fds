#!/bin/bash

# To include validation data sets to be run automatically by Validationbot,
# be sure to include the data set in ~/FDS-SMV/Validation/Process_All_Output.sh

OPTIONS="$* -y"

cd Beyler_Hood; ./Run_All.sh $OPTIONS; cd ..
cd Bittern_Sprinkler_Experiments; ./Run_All.sh $OPTIONS; cd ..
cd BRE_Spray; ./Run_All.sh $OPTIONS; cd ..
cd Bryant_Doorway; ./Run_All.sh $OPTIONS; cd ..
cd CAROLFIRE; ./Run_All.sh $OPTIONS; cd ..
cd Cup_Burner; ./Run_All.sh $OPTIONS; cd ..
cd Droplet_Evaporation; ./Run_All.sh $OPTIONS; cd ..
cd FAA_Polymers; ./Run_All.sh $OPTIONS; cd ..
cd Fleury_Heat_Flux; ./Run_All.sh $OPTIONS; cd ..
cd FM_Parallel_Panels; ./Run_All.sh $OPTIONS; cd ..
cd Hamins_Gas_Burners; ./Run_All.sh $OPTIONS; cd ..
cd Heskestad_Flame_Height; ./Run_All.sh $OPTIONS; cd ..
cd Insulation_Materials; ./Run_All.sh $OPTIONS; cd ..
cd Kashiwagi_Gasification; ./Run_All.sh $OPTIONS; cd ..
cd LEMTA_Spray; ./Run_All.sh $OPTIONS; cd ..
cd LLNL_Enclosure; ./Run_All.sh $OPTIONS; cd ..
cd NIST_He_2009; ./Run_All.sh $OPTIONS; cd ..
cd NIST_Polymers; ./Run_All.sh $OPTIONS; cd ..
cd NIST_Smoke_Alarms; ./Run_All.sh $OPTIONS; cd ..
cd NRCC_Facade; ./Run_All.sh $OPTIONS; cd ..
cd NRL_HAI; ./Run_All.sh $OPTIONS; cd ..
cd PRISME; ./Run_All.sh $OPTIONS; cd ..
cd Ranz_Marshall; ./Run_All.sh $OPTIONS; cd ..
cd Restivo_Experiment; ./Run_All.sh $OPTIONS; cd ..
cd Scaling_Pyrolysis; ./Run_All.sh $OPTIONS; cd ..
cd SP_AST; ./Run_All.sh $OPTIONS; cd ..
cd Steckler_Compartment; ./Run_All.sh $OPTIONS; cd ..
cd UL_NIST_Vents; ./Run_All.sh $OPTIONS; cd ..
cd Ulster_SBI; ./Run_All.sh $OPTIONS; cd ..
cd UMD_Polymers; ./Run_All.sh $OPTIONS; cd ..
cd UMD_Burning_Rate_Emulator; ./Run_All.sh $OPTIONS; cd ..
cd USCG_HAI; ./Run_All.sh $OPTIONS; cd ..
cd Vettori_Flat_Ceiling; ./Run_All.sh $OPTIONS; cd ..
cd Vettori_Sloped_Ceiling; ./Run_All.sh $OPTIONS; cd ..
cd VTT_Sprays; ./Run_All.sh $OPTIONS; cd ..
cd WTC; ./Run_All.sh $OPTIONS; cd ..
