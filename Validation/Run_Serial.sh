#!/bin/bash

# To include validation data sets to be run automatically by Validationbot,
# be sure to include the data set in ~/FDS-SMV/Validation/Process_All_Output.sh

OPTIONS="$* -y"

CURDIR=`pwd`

cd $CURDIR/Beyler_Hood;                   ./Run_All.sh $OPTIONS
cd $CURDIR/Bittern_Sprinkler_Experiments; ./Run_All.sh $OPTIONS
cd $CURDIR/BRE_Spray;                     ./Run_All.sh $OPTIONS
cd $CURDIR/Bryant_Doorway;                ./Run_All.sh $OPTIONS
cd $CURDIR/CAROLFIRE;                     ./Run_All.sh $OPTIONS
cd $CURDIR/CHRISTIFIRE;                   ./Run_All.sh $OPTIONS
cd $CURDIR/Cup_Burner;                    ./Run_All.sh $OPTIONS
cd $CURDIR/Droplet_Evaporation;           ./Run_All.sh $OPTIONS
cd $CURDIR/FAA_Polymers;                  ./Run_All.sh $OPTIONS
cd $CURDIR/Fleury_Heat_Flux;              ./Run_All.sh $OPTIONS
cd $CURDIR/FM_Parallel_Panels;            ./Run_All.sh $OPTIONS
cd $CURDIR/Hamins_Gas_Burners;            ./Run_All.sh $OPTIONS
cd $CURDIR/Heskestad_Flame_Height;        ./Run_All.sh $OPTIONS
cd $CURDIR/LEMTA_Spray;                   ./Run_All.sh $OPTIONS
cd $CURDIR/LLNL_Enclosure;                ./Run_All.sh $OPTIONS
cd $CURDIR/NIST_He_2009;                  ./Run_All.sh $OPTIONS
cd $CURDIR/NIST_Smoke_Alarms;             ./Run_All.sh $OPTIONS
cd $CURDIR/NIST_NRC;                      ./Run_All.sh $OPTIONS
cd $CURDIR/NRCC_Facade;                   ./Run_All.sh $OPTIONS
cd $CURDIR/NRL_HAI;                       ./Run_All.sh $OPTIONS
cd $CURDIR/PRISME;                        ./Run_All.sh $OPTIONS
cd $CURDIR/Ranz_Marshall;                 ./Run_All.sh $OPTIONS
cd $CURDIR/Restivo_Experiment;            ./Run_All.sh $OPTIONS
cd $CURDIR/SP_AST;                        ./Run_All.sh $OPTIONS
cd $CURDIR/Steckler_Compartment;          ./Run_All.sh $OPTIONS
cd $CURDIR/UL_NIST_Vents;                 ./Run_All.sh $OPTIONS
cd $CURDIR/Ulster_SBI;                    ./Run_All.sh $OPTIONS
cd $CURDIR/UMD_Polymers;                  ./Run_All.sh $OPTIONS
cd $CURDIR/UMD_Burning_Rate_Emulator;     ./Run_All.sh $OPTIONS
cd $CURDIR/USCG_HAI;                      ./Run_All.sh $OPTIONS
cd $CURDIR/Vegetation;                    ./Run_All.sh $OPTIONS
cd $CURDIR/Vettori_Flat_Ceiling;          ./Run_All.sh $OPTIONS
cd $CURDIR/Vettori_Sloped_Ceiling;        ./Run_All.sh $OPTIONS
cd $CURDIR/VTT_Sprays;                    ./Run_All.sh $OPTIONS
cd $CURDIR/WTC;                           ./Run_All.sh $OPTIONS
