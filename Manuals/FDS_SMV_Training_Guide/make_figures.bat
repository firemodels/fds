@echo off
echo creating some figures for the  Fire fighter Training guide
erase scriptfigures\*.png
cd z:\Training\Demonstrations\2Room_Ranch
smokeview -runscript ranch_01
smokeview -runscript ranch_02
smokeview -runscript ranch_03
smokeview -runscript ranch_04

cd z:\Training\MCFRS\MCFRS_Flashover
smokeview -runscript MCFRS_Flashover_00

cd z:\Training\MFRI\MFRI_Training_Tower
smokeview -runscript MFRI_Training_Tower_00
