@echo off
echo creating some figures for the  Fire fighter Training guide
erase scriptfigures\*.png
cd ..\..\Training\Demonstrations\2Room_Ranch
smokeview -runscript ranch_01
smokeview -runscript ranch_02
smokeview -runscript ranch_03
smokeview -runscript ranch_04
