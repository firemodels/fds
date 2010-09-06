@echo off
set svn_drive=d:

set BASEDIR=%CD%
set SVNROOT=%BASEDIR%\..\

set FDS=%SVNROOT%\FDS_Compilation\intel_win_32\fds5_win_32
Rem set FDS=background -u 75 -d 10 %SVNROOT%\FDS_Compilation\intel_win_32\fds5_win_32


set RUNFDS=call %SVNROOT%\Utilities\Scripts\runfds_win32.bat

echo You are about to run the Verification Test Suite.
pause > Nul

Rem -------------------------
Rem run FDS veriication cases
Rem -------------------------
echo creating FDS case list from FDS_Cases.csh
..\Utilities\Data_Processing\csh2bat FDS_Cases.csh FDS_Cases.bat

echo running FDS cases
call FDS_Cases.bat

erase FDS_Cases.bat

cd %BASEDIR%

Rem -----------------------------------------------------------------------
Rem create a text file containing the FDS5 version used to run these tests.
Rem This file is included in the smokeview user's guide
Rem -----------------------------------------------------------------------

set smvug="%CD%\..\Manuals\SMV_5_User_Guide\"
echo | fds5 2> "%smvug%\figures\fds5.version"

Rem -------------------------
Rem run FDS veriication cases
Rem -------------------------
cd %BASEDIR%\SMV_scripts
echo creating FDS case list from SMV_Cases.csh
..\..\Utilities\Data_Processing\csh2bat SMV_Cases.csh SMV_Cases.bat

call SMV_Cases.bat
erase SMV_Cases.bat

cd %BASEDIR%\SMV_scripts
call run_wui_tree_test

pause
