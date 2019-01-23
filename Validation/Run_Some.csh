#!/bin/csh -f

# To include validation data sets to be run automatically by Validationbot,
# be sure to include the data set in ~/FDS-SMV/Validation/Process_All_Output.sh

cd CSIRO_Grassland_Fires;  ./Run_All.sh -y; cd ..
cd FAA_Cargo_Compartments;  ./Run_All.sh -y; cd ..
cd NBS_Multi-Room;  ./Run_All.sh -y; cd ..
cd NRCC_Smoke_Tower;  ./Run_All.sh -y; cd ..
cd UMD_Line_Burner;  ./Run_All.sh -y; cd ..
cd UWO_Wind_Tunnel;  ./Run_All.sh -y; cd ..
