@echo off

call "%IFORT_COMPILER16%\bin\compilervars" intel64

icl -o sh2bat sh2bat.c
pause
