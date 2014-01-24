@echo off

SETLOCAL

:: ----------------------------
:: set 32 or 64 bit environment
:: ----------------------------

:: set platform=win_32
:: set platform2=ia32

set platform=win_64
set platform2=intel64

:: --------------------
:: set repository names
:: --------------------

set fdsbasename=FDS-SMV
set cfastbasename=cfast

:: --------------
:: begin smokebot
:: --------------

echo.
echo Preliminary Windows version of Smokebot
echo press any key to continue
pause>Nul

erase stage*.txt

:: -----------------
:: setup environment
:: -----------------

set CURDIR=%CD%
set OUTDIR=%CURDIR%
set svnroot=%userprofile%\%fdsbasename%
set cfastroot=%userprofile%\%cfastbasename%
set email=%svnroot%\SMV\scripts\email.bat

call "%IFORT_COMPILER14%\bin\compilervars" %platform2%
call %svnroot%\Utilities\Firebot\firebot_email_list.bat

:: -------
:: stage 0
:: -------

echo Stage 0 - Preliminaries

:: check if compilers are present

echo "" > errors.txt
ifort 1> stage0a.txt 2>&1
type stage0a.txt | find /i /c "not recognized" > count0a.txt
set /p nothaveFORTRAN=<count0a.txt
if %nothaveFORTRAN% == 1 (
  echo "***Fatal error: Fortran compiler not present"
  echo "***Fatal error: Fortran compiler not present" > errors.txt
  echo "smokebot run aborted"
  goto abort
)

icl 1> stage0b.txt 2>&1
type stage0b.txt | find /i /c "not recognized" > count0b.txt
set /p nothaveCC=<count0b.txt

:: update cfast repository

cd %cfastroot%
svn update  1> %OUTDIR%\stage0.txt 2>&1

:: build cfast

cd %cfastroot%\CFAST\intel_%platform%
erase *.obj *.mod 1> %OUTDIR%\stage0.txt 2>&1
make VPATH="../Source:../Include" INCLUDE="../Include" -f ..\makefile intel_%platform% 1>> %OUTDIR%\stage0.txt 2>&1

:: -------
:: Stage 1
:: -------

echo Stage 1 - Updating repository

cd %svnroot%
svn update 1> %OUTDIR%\stage1.txt 2>&1

:: -------
:: Stage 2
:: -------

echo Stage 2 - Building FDS (debug version)

cd %svnroot%\FDS_Compilation\intel_%platform%_db
erase *.obj *.mod 1> %OUTDIR%\stage2.txt 2>&1
make VPATH="../../FDS_Source" -f ..\makefile intel_%platform%_db 1>> %OUTDIR%\stage2.txt 2>&1

:: -------
:: Stage 3
:: -------

echo Stage 3 - Running verification cases (debug) (not implemented)

:: -------
:: Stage 4
:: -------

echo Stage 4 - Building FDS (release version)

cd %svnroot%\FDS_Compilation\intel_%platform%
erase *.obj *.mod 1> %OUTDIR%\stage4.txt 2>&1
make VPATH="../../FDS_Source" -f ..\makefile intel_%platform% 1>> %OUTDIR%\stage4.txt 2>&1

:: ----------
:: Stage 5pre
:: ----------

echo Stage 5pre - Building Smokeview utilities

echo              smokezip
cd %svnroot%\Utilities\smokezip\intel_%platform%
erase *.obj *.mod 1> %OUTDIR%\stage5pre.txt 2>&1
make -f ..\Makefile intel_%platform% 1>> %OUTDIR%\stage5pre.txt 2>&1

echo              smokediff
cd %svnroot%\Utilities\smokediff\intel_%platform%
erase *.obj *.mod 1>> %OUTDIR%\stage5pre.txt 2>&1
make -f ..\Makefile intel_%platform% 1>> %OUTDIR%\stage5pre.txt 2>&1

echo              wind2fds
cd %svnroot%\Utilities\wind2fds\intel_%platform%
erase *.obj *.mod 1>> %OUTDIR%\stage5pre.txt 2>&1
make -f ..\Makefile intel_%platform% 1>> %OUTDIR%\stage5pre.txt 2>&1

:: -------
:: Stage 5
:: -------

echo Stage 5 - Running verification cases (release)

cd %svnroot%\Verification\scripts
call run_smv_cases 1> %OUTDIR%\stage5.txt 2>&1

:: --------
:: Stage 6a
:: --------

echo Stage 6a - Building Smokeview (debug version)

cd %svnroot%\SMV\Build\intel_%platform%_db
erase *.obj *.mod 1> %OUTDIR%\stage6a.txt 2>&1
make -f ..\Makefile intel_%platform%_db 1>> %OUTDIR%\stage6a.txt 2>&1

:: --------
:: Stage 6b
:: --------

echo Stage 6b - Making Smokeview pictures (debug mode) (not implemented)

:: --------
:: Stage 6c
:: --------

echo Stage 6c - Building Smokeview (release version)

cd %svnroot%\SMV\Build\intel_%platform%
erase *.obj *.mod 1> %OUTDIR%\stage6c.txt 2>&1
make -f ..\Makefile intel_%platform% 1>> %OUTDIR%\stage6c.txt 2>&1

:: --------
:: Stage 6d
:: --------

echo Stage 6d - Making Smokeview pictures (release mode)

cd %svnroot%\Verification\scripts
:: call MAKE_SMV_pictures 1> %OUTDIR%\stage6d.txt 2>&1
call MAKE_SMV_pictures

:: -------
:: Stage 8
:: -------

echo Stage 8 - Building Smokeview guides

cd %svnroot%\Manuals\SMV_User_Guide
echo                 User
call make_guide 1> %OUTDIR%\stage8.txt 2>&1

cd %svnroot%\Manuals\SMV_Technical_Reference_Guide
echo                 Technical Reference
call make_guide 1>> %OUTDIR%\stage8.txt 2>&1

cd %svnroot%\Manuals\SMV_Verification_Guide
echo                 Verification
call make_guide 1>> %OUTDIR%\stage8.txt 2>&1

echo smokebot build success on %COMPUTERNAME% > info.txt

cd %CURDIR%

:: -------
:: Wrap up
:: -------

%email% %mailToSMV% "smokebot build success on %COMPUTERNAME%" info.txt

echo smokebot_win completed
pause
goto eof

:abort
%email% %mailToSMV% "smokebot build failure on %COMPUTERNAME%" errors.txt



:eof