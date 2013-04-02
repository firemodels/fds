@echo off

call "%IFORT_COMPILER13%\bin\compilervars" intel64 vs2008

icl -o sh2bat sh2bat.c
pause
