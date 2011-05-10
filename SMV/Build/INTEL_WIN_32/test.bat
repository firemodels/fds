@echo off

echo ***************** setting up fortran environment ********************
call "%IFORT_COMPILER11%\bin\ifortvars" ia32 vs2008
pause
echo ***************** setting up C environment ***************************
call "%IFORT_COMPILER11%\bin\iclvars" ia32 vs2008
pause
echo ***************** testing fortran ***************************
ifort
pause
echo ***************** testing C **************************
icl
pause
echo ***************** testing make ***********************
make
pause