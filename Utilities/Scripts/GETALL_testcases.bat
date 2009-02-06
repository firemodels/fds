@echo off

Rem get test cases for Windows and Linux/OSX 

set revision=3222
set REPOS=d:\fds-smv
set LREPOS=FDS-SMV
set logon=gforney@acrux.cfr.nist.gov

Rem -------- should not need to edit below -----------

set scriptdir=%LREPOS%/Utilities/Scripts
set testdir=FDS_Test_cases_%revision%

Rem get test cases for Linux/OSX systems

plink %logon% %scriptdir%/GET_testcases.csh %revision%
pscp %logon%:%scriptdir%/to_google/%testdir%.tar.gz to_google\.


Rem get test cases for Windows systems

call GET_testcases %revision%

pause
