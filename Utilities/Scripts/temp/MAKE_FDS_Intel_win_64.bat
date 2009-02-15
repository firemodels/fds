@echo off

Rem Batch file used to build a 64 bit version of FDS

Rem location of batch files used to set up Intel compilation environment
set intelbin=c:\bin

Rem --------- Should not need to edit below -------------

call %intelbin%\iclvars intel64
call %intelbin%\ifortvars intel64

cd ..\Makefile\Intel_Win_64

Rem If you want to do a full compile then remove following two Rem's
Rem erase *.obj
Rem erase *.mod
%intelbin%\make VPATH="../../../FDS_Source" -f ..\makefile intel_win_64

cd ..\..\Scripts


pause