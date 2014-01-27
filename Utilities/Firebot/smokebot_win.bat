@echo off

SETLOCAL

:: ----------------------------
:: set 32 or 64 bit environment
:: ----------------------------

:: set size=32
:: set compile_platform=ia32

set size=64
set compile_platform=intel64

:: --------------------
:: set repository names
:: --------------------

set fdsbasename=FDS-SMV
set cfastbasename=cfast

:: --------------
:: begin smokebot
:: --------------

erase stage*.txt

:: -----------------
:: setup environment
:: -----------------

set CURDIR=%CD%
set OUTDIR=%CURDIR%
set svnroot=%userprofile%\%fdsbasename%
set cfastroot=%userprofile%\%cfastbasename%
set email=%svnroot%\SMV\scripts\email.bat

set mpichinc="c:\mpich\mpich2_%size%\include"
set mpichlib="c:\mpich\mpich2_%size%\lib\fmpich2.lib"

set errorfile=%OUTDIR%\errors.txt
set warningfile=%OUTDIR%\warnings.txt
set infofile=%OUTDIR%\info.txt
set revisionfile=%OUTDIR%\revision.txt

set havewarnings=0
set haveCC=1

set emailexe=%userprofile%\bin\mailsend.exe

call "%IFORT_COMPILER14%\bin\compilervars" %compile_platform% 1> Nul 2>&1
call %svnroot%\Utilities\Firebot\firebot_email_list.bat

:: -------
:: stage 0
:: -------

echo Stage 0 - Preliminaries

:: check if compilers are present

echo. > %errorfile%
echo. > %warningfile%

ifort 1> stage0a.txt 2>&1
type stage0a.txt | find /i /c "not recognized" > count0a.txt
set /p nothaveFORTRAN=<count0a.txt
if %nothaveFORTRAN% == 1 (
  echo "***Fatal error: Fortran compiler not present"
  echo "***Fatal error: Fortran compiler not present" > %errorfile%
  echo "smokebot run aborted"
  call :output_abort_message
  exit /b 1
)
echo             found Fortran

icl 1> stage0b.txt 2>&1
type stage0b.txt | find /i /c "not recognized" > count0b.txt
set /p nothaveCC=<count0b.txt
if %nothaveCC% == 1 (
  set haveCC=0
  echo "***Warning: C/C++ compiler not found - using installed Smokeview to generate images"
) else (
  echo             found C/C++
)

if NOT exist %emailexe% (
  echo ***warning: email client not found.   
  echo             Smokebot messages will only be sent to the console.
) else (
  echo             found mailsend
)

call :is_file_installed pdflatex|| exit /b 1
echo             found pdflatex

:: update cfast repository

echo             updating cfast repository
cd %cfastroot%
svn update  1> %OUTDIR%\stage0.txt 2>&1

:: build cfast

echo             building cfast
cd %cfastroot%\CFAST\intel_win_%size%
erase *.obj *.mod *.exe 1> %OUTDIR%\stage0.txt 2>&1
make VPATH="../Source:../Include" INCLUDE="../Include" -f ..\makefile intel_win_%size% 1>> %OUTDIR%\stage0.txt 2>&1
call :does_file_exist cfast6_win_%size%.exe %OUTDIR%\stage0.txt|| exit /b 1

:: -------
:: Stage 1
:: -------

echo Stage 1 - Updating FDS/Smokeview repository

cd %svnroot%
svn update 1> %OUTDIR%\stage1.txt 2>&1

svn info | find /i "Revision" > %revisionfile%
set /p revision=<%revisionfile%

:: -------
:: Stage 2
:: -------

echo Stage 2 - Building FDS (debug version)
echo             serial
cd %svnroot%\FDS_Compilation\intel_win_%size%_db
erase *.obj *.mod *.exe 1> %OUTDIR%\stage2a.txt 2>&1
make -j4 VPATH="../../FDS_Source" -f ..\makefile intel_win_%size%_db 1>> %OUTDIR%\stage2a.txt 2>&1

call :does_file_exist fds_win_%size%_db.exe %OUTDIR%\stage2a.txt|| exit /b 1
call :find_fds_string "warning" %OUTDIR%\stage2a.txt

echo             parallel
cd %svnroot%\FDS_Compilation\mpi_intel_win_%size%_db
erase *.obj *.mod *.exe 1> %OUTDIR%\stage2b.txt 2>&1
make -j4 MPIINCLUDE=%mpichinc% MPILIB=%mpichlib% VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_%size%_db 1>> %OUTDIR%\stage2b.txt 2>&1

call :does_file_exist fds_mpi_win_%size%_db.exe %OUTDIR%\stage2b.txt|| exit /b 1
call :find_fds_string "warning" %OUTDIR%\stage2b.txt

:: -------
:: Stage 3
:: -------

echo Stage 3 - Running verification cases (debug) (not implemented)

:: -------
:: Stage 4
:: -------

echo Stage 4 - Building FDS (release version)

echo             serial
cd %svnroot%\FDS_Compilation\intel_win_%size%
erase *.obj *.mod *.exe 1> %OUTDIR%\stage4a.txt 2>&1
make -j4 VPATH="../../FDS_Source" -f ..\makefile intel_win_%size% 1>> %OUTDIR%\stage4a.txt 2>&1

call :does_file_exist fds_win_%size%.exe %OUTDIR%\stage4a.txt|| exit /b 1
call :find_fds_string "warning" %OUTDIR%\stage4a.txt

echo             parallel
cd %svnroot%\FDS_Compilation\mpi_intel_win_%size%
erase *.obj *.mod *.exe 1> %OUTDIR%\stage4b.txt 2>&1
make -j4 MPIINCLUDE=%mpichinc% MPILIB=%mpichlib% VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_%size%  1>> %OUTDIR%\stage4b.txt 2>&1

call :does_file_exist fds_mpi_win_%size%.exe %OUTDIR%\stage4b.txt|| exit /b 1
call :find_fds_string "warning" %OUTDIR%\stage4b.txt

:: ----------
:: Stage 5pre
:: ----------

echo Stage 5pre - Building Smokeview utilities

echo                fds2ascii
cd %svnroot%\Utilities\fds2ascii\intel_win_%size%
erase *.obj *.mod *.exe 1> %OUTDIR%\stage5pre.txt 2>&1
make -f ..\Makefile intel_win_%size% 1>> %OUTDIR%\stage5pre.txt 2>&1
call :does_file_exist fds2ascii_win_%size%.exe %OUTDIR%\stage5pre.txt|| exit /b 1

if %haveCC% == 1 (
  echo                smokediff
  cd %svnroot%\Utilities\smokediff\intel_win_%size%
  erase *.obj *.mod *.exe 1>> %OUTDIR%\stage5pre.txt 2>&1
  make -f ..\Makefile intel_win_%size% 1>> %OUTDIR%\stage5pre.txt 2>&1
  call :does_file_exist smokediff_win_%size%.exe %OUTDIR%\stage5pre.txt

  echo                smokezip
  cd %svnroot%\Utilities\smokezip\intel_win_%size%
  erase *.obj *.mod *.exe 1>> %OUTDIR%\stage5pre.txt 2>&1
  make -f ..\Makefile intel_win_%size% 1>> %OUTDIR%\stage5pre.txt 2>&1
  call :does_file_exist smokezip_win_%size%.exe %OUTDIR%\stage5pre.txt|| exit /b 1

  echo                wind2fds
  cd %svnroot%\Utilities\wind2fds\intel_win_%size%
  erase *.obj *.mod *.exe 1>> %OUTDIR%\stage5pre.txt 2>&1
  make -f ..\Makefile intel_win_%size% 1>> %OUTDIR%\stage5pre.txt 2>&1
  call :does_file_exist wind2fds_win_%size%.exe %OUTDIR%\stage5pre.txt|| exit /b 1
) else (
  call :is_file_installed smokediff|| exit /b 1
  echo                smokediff not built, is installed
  call :is_file_installed smokezip|| exit /b 1
  echo                smokezip not built, is installed
  call :is_file_installed wind2fds|| exit /b 1
  echo                wind2fds not built, is installed
)

:: -------
:: Stage 5
:: -------

echo Stage 5 - Running verification cases (release)

cd %svnroot%\Verification\scripts
call Run_SMV_cases %size% 1> %OUTDIR%\stage5.txt 2>&1

:: --------
:: Stage 6a
:: --------

echo Stage 6a - Building Smokeview (debug version)

cd %svnroot%\SMV\Build\intel_win_%size%_db
erase *.obj *.mod *.exe 1> %OUTDIR%\stage6a.txt 2>&1
make -f ..\Makefile intel_win_%size%_db 1>> %OUTDIR%\stage6a.txt 2>&1

call :does_file_exist smokeview_win_%size%_db.exe %OUTDIR%\stage6a.txt|| exit /b 1
call :find_string "warning" %OUTDIR%\stage6a.txt

:: --------
:: Stage 6b
:: --------

echo Stage 6b - Making Smokeview pictures (debug mode) (not implemented)

:: --------
:: Stage 6c
:: --------

echo Stage 6c - Building Smokeview (release version)

cd %svnroot%\SMV\Build\intel_win_%size%
erase *.obj *.mod *.exe 1> %OUTDIR%\stage6c.txt 2>&1
make -f ..\Makefile intel_win_%size% 1>> %OUTDIR%\stage6c.txt 2>&1

call :does_file_exist smokeview_win_%size%.exe %OUTDIR%\stage6c.txt|| exit /b 1
call :find_string "warning" %OUTDIR%\stage6c.txt

:: --------
:: Stage 6d
:: --------

echo Stage 6d - Making Smokeview pictures (release mode)

cd %svnroot%\Verification\scripts
call MAKE_SMV_pictures %size% 1> %OUTDIR%\stage6d.txt 2>&1

call :find_string "error" %OUTDIR%\stage6d.txt

:: -------
:: Stage 8
:: -------

echo Stage 8 - Building Smokeview guides

cd %svnroot%\Manuals\SMV_User_Guide
echo             User
call make_guide 1> %OUTDIR%\stage8.txt 2>&1

cd %svnroot%\Manuals\SMV_Technical_Reference_Guide
echo             Technical Reference
call make_guide 1>> %OUTDIR%\stage8.txt 2>&1

cd %svnroot%\Manuals\SMV_Verification_Guide
echo             Verification
call make_guide 1>> %OUTDIR%\stage8.txt 2>&1

echo smokebot build success on %COMPUTERNAME% > %infofile%

cd %CURDIR%

:: -------
:: Wrap up
:: -------

if exist %emailexe% (
  if %havewarnings% == 0 (
    call %email% %mailToSMV% "smokebot build success on %COMPUTERNAME%! %revision%" %infofile%  
  )
  if %havewarnings% GTR 0 (
    %email% %mailToSMV% "smokebot build success with warnings on %COMPUTERNAME% %revision%" %warningfile%
  )
)

echo smokebot_win completed
goto eof

:output_abort_message
  echo "***Fatal error: smokebot build failure on %COMPUTERNAME% %revision%"
  if exist %email% (
    call %email% %mailToSMV% "smokebot build failure on %COMPUTERNAME% %revision%" %errorfile%
  )
exit /b

:: -----------------------------------------
:is_file_installed
:: -----------------------------------------

  set program=%1
  %program% -help 1>> exist.txt 2>&1
  type exist.txt | find /i /c "not recognized" > count.txt
  set /p nothave=<count.txt
  if %nothave% == 1 (
    echo "***Fatal error: %program% not present"
    echo "***Fatal error: %program% not present" > %errorfile%
    echo "smokebot run aborted"
    call :output_abort_message
    exit /b 1
  )
  exit /b 0

:: -----------------------------------------
  :does_file_exist
:: -----------------------------------------

set file=%1
set outputfile=%2

if NOT exist %file% (
  echo ***fatal error: problem building %file%. Aborting smokebot
  type %outputfile% >> %errorfile%
  call :output_abort_message
  exit /b 1
)
exit /b 0

:: -----------------------------------------
  :find_string
:: -----------------------------------------

set search_string=%1
set search_file=%2

findstr /I %search_string% %search_file% > %OUTDIR%\warning.txt
type %OUTDIR%\warning.txt | find /c ":"> %OUTDIR%\nwarning.txt
set /p nwarnings=<%OUTDIR%\nwarning.txt
if %nwarnings% GTR 0 (
  type %OUTDIR%\warning.txt >> %warningfile%
  set havewarnings=1
)
exit /b

:: -----------------------------------------
  :find_fds_string
:: -----------------------------------------

set search_string=%1
set search_file=%2

findstr /I %search_string% %search_file% | find /V "mpif.h" | find /V "commands for target" > %OUTDIR%\warning.txt
type %OUTDIR%\warning.txt | find /c ":"> %OUTDIR%\nwarning.txt
set /p nwarnings=<%OUTDIR%\nwarning.txt
if %nwarnings% GTR 0 (
  type %OUTDIR%\warning.txt >> %warningfile%
  set havewarnings=1
)
exit /b


:eof
cd %CURDIR%
pause
exit
