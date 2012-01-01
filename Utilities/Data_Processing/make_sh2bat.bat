@echo off

call "%IFORT_COMPILER12%\bin\iclvars" intel64 vs2008
call "%IFORT_COMPILER12%\bin\iclvars" ia32 vs2008

icl -o sh2bat sh2bat.c
pause
