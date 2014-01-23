@echo off

SETLOCAL
set platform=win_64
set fdsbasename=FDS-SMV
set cfastbasename=cfast

echo.
echo Preliminary Windows version of Smokebot
echo press any key to continue
pause>Nul

Rem ---------------------------------------------
Rem set up environment variables used by smokebot
Rem ---------------------------------------------

set CURDIR=%CD%
erase stage*.txt

set OUTDIR=%CURDIR%

set svnroot=%userprofile%\%fdsbasename%
set cfastroot=%userprofile%\%cfastbasename%

Rem ---------------------------------------
echo Stage 0 - Building cfast
Rem ---------------------------------------

cd %cfastroot%\CFAST\intel_%platform%
erase *.obj *.mod 1> %OUTDIR%\stage0.txt 2>&1
call make_cfast 1>> %OUTDIR%\stage0.txt 2>&1

Rem -------------------------
echo Stage 1 - Updating repository
Rem -------------------------

cd %svnroot%
svn update 1> %OUTDIR%\stage1.txt 2>&1

Rem --------------------------
echo Stage 2 - Building FDS (debug version)
Rem --------------------------

cd %svnroot%\FDS_Compilation\intel_%platform%_db
erase *.obj *.mod 1> %OUTDIR%\stage2.txt 2>&1
call make_fds.bat 1>> %OUTDIR%\stage2.txt 2>&1

Rem ----------------------------------------
echo Stage 3 - Running verification cases (debug) (not implemented)
Rem ----------------------------------------

Rem ---------------------------  
echo Stage 4 - Building FDS (release version)
Rem ---------------------------  

cd %svnroot%\FDS_Compilation\intel_%platform%
erase *.obj *.mod 1> %OUTDIR%\stage4.txt 2>&1
call make_fds.bat 1>> %OUTDIR%\stage4.txt 2>&1

Rem -----------------------------------  
echo Stage 5pre - Building Smokeview utilities
Rem -----------------------------------

echo              smokezip
cd %svnroot%\Utilities\smokezip\intel_%platform%
call make_zip.bat 1> %OUTDIR%\stage5pre.txt 2>&1

echo              smokediff
cd %svnroot%\Utilities\smokediff\intel_%platform%
call make_diff.bat 1>> %OUTDIR%\stage5pre.txt 2>&1

echo              wind2fds
cd %svnroot%\Utilities\wind2fds\intel_%platform%
call make_wind.bat 1>> %OUTDIR%\stage5pre.txt 2>&1

Rem ------------------------------------------
echo Stage 5 - Running verification cases (release)
Rem ------------------------------------------

cd %svnroot%\Verification\scripts
call run_smv_cases 1> %OUTDIR%\stage5.txt 2>&1

Rem ---------------------------------
echo Stage 6a - Building Smokeview (debug version)
Rem ---------------------------------

cd %svnroot%\SMV\Build\intel_%platform%_db
call make_smv.bat 1> %OUTDIR%\stage6a.txt 2>&1

Rem -----------------------------------------
echo Stage 6b - Making Smokeview pictures (debug mode) (not implemented)
Rem -----------------------------------------

Rem ---------------------------------
echo Stage 6c - Building Smokeview (release version)
Rem ---------------------------------

cd %svnroot%\SMV\Build\intel_%platform%
call make_smv.bat 1> %OUTDIR%\stage6c.txt 2>&1

Rem -----------------------------------------
echo Stage 6d - Making Smokeview pictures (release mode)
Rem -----------------------------------------

cd %svnroot%\Verification\scripts
Rem call MAKE_SMV_pictures 1> %OUTDIR%\stage6d.txt 2>&1
call MAKE_SMV_pictures

Rem ---------------------------------
echo Stage 8 - Building Smokeview guides
Rem ---------------------------------

cd %svnroot%\Manuals\SMV_User_Guide
echo                 User
call make_guide 1> %OUTDIR%\stage8.txt 2>&1

cd %svnroot%\Manuals\SMV_Technical_Reference_Guide
echo                 Technical Reference
call make_guide 1>> %OUTDIR%\stage8.txt 2>&1

cd %svnroot%\Manuals\SMV_Verification_Guide
echo                 Verification
call make_guide 1>> %OUTDIR%\stage8.txt 2>&1

cd %CURDIR%

Rem -------
Rem Wrap up
Rem -------

echo smokebot_win completed
pause