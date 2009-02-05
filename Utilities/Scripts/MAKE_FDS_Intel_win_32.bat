@echo off

Rem Batch file used to build a 32 bit version of FDS

Rem location of batch files used to set up Intel compilation environment
set intelbin=d:\bin

Rem --------- Should not need to edit below -------------

call %intelbin%\iclvars ia32
call %intelbin%\ifortvars ia32

cd ..\Makefile\Intel_Win_32

Rem If you want to do a full compile then remove following two Rem's
Rem erase *.obj
Rem erase *.mod
make VPATH="../../../FDS_Source" -f ..\makefile intel_win_32

cd ..\..\Scripts


pause