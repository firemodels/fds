#!/bin/csh -f
echo "****** Copy_csv *******"
./Copy_csv.csh
echo "****** Make_ascii_csv *******"
./Make_ascii_csv.csh
echo "****** Merge_csv *******"
./Merge_csv.csh
echo "****** Merge_door_csv *******"
./Merge_door_csv.csh
echo "****** Merge_SL_csv *******"
./Merge_SL_csv.csh
