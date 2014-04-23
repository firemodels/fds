@echo off
set reduced=%1
if [%reduced%] == [] (
  set reduced=0
)


:: -------------------------------------------------------------
::                         set 32 or 64 bit environment
:: -------------------------------------------------------------

:: set size=32
:: set compile_platform=ia32

set size=64
set compile_platform=intel64

set OPENMP=
:: to use non-openmp fds comment the following two lines
set OPENMP=openmp_
set OMP_NUM_THREADS=2

:: -------------------------------------------------------------
::                         set repository names
:: -------------------------------------------------------------

set fdsbasename=FDS-SMV
set cfastbasename=cfast

:: -------------------------------------------------------------
::                         setup environment
:: -------------------------------------------------------------

set CURDIR=%CD%

if not exist output mkdir output

set OUTDIR=%CURDIR%\output

erase %OUTDIR%\*.txt 1> Nul 2>&1

set svnroot=%userprofile%\%fdsbasename%
set cfastroot=%userprofile%\%cfastbasename%
set email=%svnroot%\SMV\scripts\email.bat

set mpichinc="c:\mpich\mpich2_%size%\include"
set mpichlib="c:\mpich\mpich2_%size%\lib\fmpich2.lib"

set errorlog=%OUTDIR%\stage_errors.txt
set warninglog=%OUTDIR%\stage_warnings.txt
set errorwarninglog=%OUTDIR%\stage_errorswarnings.txt
set infofile=%OUTDIR%\stage_info.txt
set revisionfile=%OUTDIR%\revision.txt
set stagestatus=%OUTDIR%\stage_status.log

set fromsummarydir=%svnroot%\Manuals\SMV_Summary
set tosummarydir="%SMOKEBOT_SUMMARY_DIR%"

set bundle_win32=%svnroot%\Utilities\Scripts\BUNDLE_win32.bat
set bundle_win64=%svnroot%\Utilities\Scripts\BUNDLE_win64.bat
set bundle_linux32=%svnroot%\Utilities\Scripts\BUNDLE_linux32.bat
set bundle_linux64=%svnroot%\Utilities\Scripts\BUNDLE_linux64.bat
set bundle_osx32=%svnroot%\Utilities\Scripts\BUNDLE_osx32.bat
set bundle_osx64=%svnroot%\Utilities\Scripts\BUNDLE_osx64.bat

set haveerrors=0
set havewarnings=0
set haveCC=1

set emailexe=%userprofile%\bin\mailsend.exe
set gettimeexe=%svnroot%\Utilities\get_time\intel_win_64\get_time.exe

date /t > %OUTDIR%\starttime.txt
set /p startdate=<%OUTDIR%\starttime.txt
time /t > %OUTDIR%\starttime.txt
set /p starttime=<%OUTDIR%\starttime.txt

call "%IFORT_COMPILER14%\bin\compilervars" %compile_platform% 1> Nul 2>&1
call %svnroot%\Utilities\Firebot\firebot_email_list.bat

:: -------------------------------------------------------------
::                           stage 0
:: -------------------------------------------------------------

echo Stage 0 - Preliminaries

:: check if compilers are present

echo. > %errorlog%
echo. > %warninglog%
echo. > %stagestatus%

call :is_file_installed %gettimeexe%|| exit /b 1
echo             found get_time

call :GET_TIME
set TIME_beg=%current_time%

call :GET_TIME
set PRELIM_beg=%current_time% 

ifort 1> %OUTDIR%\stage0a.txt 2>&1
type %OUTDIR%\stage0a.txt | find /i /c "not recognized" > %OUTDIR%\stage_count0a.txt
set /p nothaveFORTRAN=<%OUTDIR%\stage_count0a.txt
if %nothaveFORTRAN% == 1 (
  echo "***Fatal error: Fortran compiler not present"
  echo "***Fatal error: Fortran compiler not present" > %errorlog%
  echo "smokebot run aborted"
  call :output_abort_message
  exit /b 1
)
echo             found Fortran

icl 1> %OUTDIR%\stage0b.txt 2>&1
type %OUTDIR%\stage0b.txt | find /i /c "not recognized" > %OUTDIR%\stage_count0b.txt
set /p nothaveCC=<%OUTDIR%\stage_count0b.txt
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

:: update FDS/Smokeview repository

echo             updating FDS/Smokeview repository

cd %svnroot%
svn update 1>> %OUTDIR%\stage0.txt 2>&1

svn info | find /i "Revision" > %revisionfile%
set /p revision=<%revisionfile%

:: build cfast

echo             building cfast
cd %cfastroot%\CFAST\intel_win_%size%
erase *.obj *.mod *.exe 1>> %OUTDIR%\stage0.txt 2>&1
make VPATH="../Source:../Include" INCLUDE="../Include" -f ..\makefile intel_win_%size% 1>> %OUTDIR%\stage0.txt 2>&1
call :does_file_exist cfast6_win_%size%.exe %OUTDIR%\stage0.txt|| exit /b 1

echo             building geomtest
cd %svnroot%\SMV\source\geomtest\intel_win_%size%
erase *.obj *.mod *.exe 1>> %OUTDIR%\stage0.txt 2>&1
make VPATH=".." -f ..\makefile intel_win_%size% 1>> %OUTDIR%\stage0.txt 2>&1
call :does_file_exist geomtest.exe %OUTDIR%\stage0.txt|| exit /b 1

call :GET_TIME
set PRELIM_end=%current_time%
call :GET_DURATION PRELIM %PRELIM_beg% %PRELIM_end%
set DIFF_PRELIM=%duration%

:: -------------------------------------------------------------
::                           stage 1
:: -------------------------------------------------------------

call :GET_TIME
set BUILDFDS_beg=%current_time% 
echo Stage 1 - Building FDS
if %reduced% == 1 goto skip_fds_debug

echo             serial debug

cd %svnroot%\FDS_Compilation\%OPENMP%intel_win_%size%_db
erase *.obj *.mod *.exe 1> %OUTDIR%\stage1a.txt 2>&1
make VPATH="../../FDS_Source" -f ..\makefile %OPENMP%intel_win_%size%_db 1>> %OUTDIR%\stage1a.txt 2>&1

call :does_file_exist fds_%OPENMP%win_%size%_db.exe %OUTDIR%\stage1a.txt|| exit /b 1
call :find_fds_warnings "warning" %OUTDIR%\stage1a.txt "Stage 1a"

echo             parallel debug

cd %svnroot%\FDS_Compilation\mpi_intel_win_%size%_db
erase *.obj *.mod *.exe 1> %OUTDIR%\stage1b.txt 2>&1
make MPIINCLUDE=%mpichinc% MPILIB=%mpichlib% VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_%size%_db 1>> %OUTDIR%\stage1b.txt 2>&1

call :does_file_exist fds_mpi_win_%size%_db.exe %OUTDIR%\stage1b.txt|| exit /b 1
call :find_fds_warnings "warning" %OUTDIR%\stage1b.txt "Stage 1b"

:skip_fds_debug

echo             serial release

cd %svnroot%\FDS_Compilation\%OPENMP%intel_win_%size%
erase *.obj *.mod *.exe 1> %OUTDIR%\stage1c.txt 2>&1
make VPATH="../../FDS_Source" -f ..\makefile %OPENMP%intel_win_%size% 1>> %OUTDIR%\stage1c.txt 2>&1

call :does_file_exist fds_%OPENMP%win_%size%.exe %OUTDIR%\stage1c.txt|| exit /b 1
call :find_fds_warnings "warning" %OUTDIR%\stage1c.txt "Stage 1c"

if %reduced% == 1 goto skip_fds_parallel

echo             parallel release

cd %svnroot%\FDS_Compilation\mpi_intel_win_%size%
erase *.obj *.mod *.exe 1> %OUTDIR%\stage1d.txt 2>&1
make MPIINCLUDE=%mpichinc% MPILIB=%mpichlib% VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_%size%  1>> %OUTDIR%\stage1d.txt 2>&1

call :does_file_exist fds_mpi_win_%size%.exe %OUTDIR%\stage1d.txt|| exit /b 1
call :find_fds_warnings "warning" %OUTDIR%\stage1d.txt "Stage 1d"

:skip_fds_parallel

call :GET_TIME
set BUILDFDS_end=%current_time%
call :GET_DURATION BUILDFDS %BUILDFDS_beg% %BUILDFDS_end%
set DIFF_BUILDFDS=%duration%

:: -------------------------------------------------------------
::                           stage 2
:: -------------------------------------------------------------

call :GET_TIME
set BUILDSMVUTIL_beg=%current_time% 
echo Stage 2 - Building Smokeview

if %reduced% == 1 goto skip_smokeview_debug

echo             debug

cd %svnroot%\SMV\Build\intel_win_%size%
erase *.obj *.mod *.exe smokeview_win_%size%_db.exe 1> %OUTDIR%\stage2a.txt 2>&1
make -f ..\Makefile intel_win_%size%_db 1>> %OUTDIR%\stage2a.txt 2>&1

call :does_file_exist smokeview_win_%size%_db.exe %OUTDIR%\stage2a.txt|| exit /b 1
call :find_smokeview_warnings "warning" %OUTDIR%\stage2a.txt "Stage 2a"

:skip_smokeview_debug

echo             release

cd %svnroot%\SMV\Build\intel_win_%size%
erase *.obj *.mod smokeview_win_%size%.exe 1> %OUTDIR%\stage2b.txt 2>&1
make -f ..\Makefile intel_win_%size% 1>> %OUTDIR%\stage2b.txt 2>&1

call :does_file_exist smokeview_win_%size%.exe %OUTDIR%\stage2b.txt|| aexit /b 1
call :find_smokeview_warnings "warning" %OUTDIR%\stage2b.txt "Stage 2b"

:: -------------------------------------------------------------
::                           stage 3
:: -------------------------------------------------------------

echo Stage 3 - Building FDS/Smokeview utilities

if %reduced% == 1 goto skip_fds2ascii

echo             fds2ascii
cd %svnroot%\Utilities\fds2ascii\intel_win_%size%
erase *.obj *.mod *.exe 1> %OUTDIR%\stage3c.txt 2>&1
make -f ..\Makefile intel_win_%size% 1>> %OUTDIR%\stage3.txt 2>&1
call :does_file_exist fds2ascii_win_%size%.exe %OUTDIR%\stage3.txt|| exit /b 1

:skip_fds2ascii

if %haveCC% == 1 (
  echo             smokediff
  cd %svnroot%\Utilities\smokediff\intel_win_%size%
  erase *.obj *.mod *.exe 1>> %OUTDIR%\stage3.txt 2>&1
  make -f ..\Makefile intel_win_%size% 1>> %OUTDIR%\stage3.txt 2>&1
  call :does_file_exist smokediff_win_%size%.exe %OUTDIR%\stage3.txt

  echo             smokezip
  cd %svnroot%\Utilities\smokezip\intel_win_%size%
  erase *.obj *.mod *.exe 1>> %OUTDIR%\stage3.txt 2>&1
  make -f ..\Makefile intel_win_%size% 1>> %OUTDIR%\stage3.txt 2>&1
  call :does_file_exist smokezip_win_%size%.exe %OUTDIR%\stage3.txt|| exit /b 1

  echo             wind2fds
  cd %svnroot%\Utilities\wind2fds\intel_win_%size%
  erase *.obj *.mod *.exe 1>> %OUTDIR%\stage3.txt 2>&1
  make -f ..\Makefile intel_win_%size% 1>> %OUTDIR%\stage3.txt 2>&1
  call :does_file_exist wind2fds_win_%size%.exe %OUTDIR%\stage3.txt|| exit /b 1
) else (
  call :is_file_installed smokediff|| exit /b 1
  echo             smokediff not built, using installed version
  call :is_file_installed smokezip|| exit /b 1
  echo             smokezip not built, using installed version
  call :is_file_installed wind2fds|| exit /b 1
  echo             wind2fds not built, using installed version
)

call :GET_TIME
set BUILDSMVUTIL_end=%current_time% 
call :GET_DURATION BUILDSMVUTIL %BUILDSMVUTIL_beg% %BUILDSMVUTIL_end%
set DIFF_BUILDSMVUTIL=%duration%

:: -------------------------------------------------------------
::                           stage 4
:: -------------------------------------------------------------

call :GET_TIME
set RUNVV_beg=%current_time% 

echo Stage 4 - Running verification cases

cd %svnroot%\Verification\scripts
call Run_SMV_cases %size% 1> %OUTDIR%\stage4.txt 2>&1

call :find_smokeview_warnings "error" %OUTDIR%\stage4.txt "Stage 4"

call :GET_TIME
set RUNVV_end=%current_time% 
call :GET_DURATION RUNVV %RUNVV_beg% %RUNVV_end%
set DIFF_RUNVV=%duration%

:: -------------------------------------------------------------
::                           stage 5
:: -------------------------------------------------------------

call :GET_TIME
set MAKEPICS_beg=%current_time% 
echo Stage 5 - Making Smokeview pictures

cd %svnroot%\Verification\scripts
call MAKE_SMV_pictures %size% 1> %OUTDIR%\stage5.txt 2>&1

call :find_smokeview_warnings "error" %OUTDIR%\stage5.txt "Stage 5"

call :GET_TIME
set MAKEPICS_end=%current_time% 
call :GET_DURATION MAKEPICS %MAKEPICS_beg% %MAKEPICS_end%
set DIFF_MAKEPICS=%duration%

:: -------------------------------------------------------------
::                           stage 6
:: -------------------------------------------------------------

call :GET_TIME
set MAKEGUIDES_beg=%current_time% 
echo Stage 6 - Building Smokeview guides

echo             Technical Reference
call :build_guide SMV_Technical_Reference_Guide %svnroot%\Manuals\SMV_Technical_Reference_Guide 1>> %OUTDIR%\stage6.txt 2>&1

echo             Verification
call :build_guide SMV_Verification_Guide %svnroot%\Manuals\SMV_Verification_Guide 1>> %OUTDIR%\stage6.txt 2>&1

echo             User
call :build_guide SMV_User_Guide %svnroot%\Manuals\SMV_User_Guide 1>> %OUTDIR%\stage6.txt 2>&1

call :GET_TIME
set MAKEGUIDES_end=%current_time%
call :GET_DURATION MAKEGUIDES %MAKEGUIDES_beg% %MAKEGUIDES_end%
set DIFF_MAKEGUIDES=%duration%

call :GET_TIME
set TIME_end=%current_time% 
call :GET_DURATION TOTALTIME %TIME_beg% %TIME_end%
set DIFF_TIME=%duration%

:: -------------------------------------------------------------
::                           wrap up
:: -------------------------------------------------------------

date /t > %OUTDIR%\stoptime.txt
set /p stopdate=<%OUTDIR%\stoptime.txt
time /t > %OUTDIR%\stoptime.txt
set /p stoptime=<%OUTDIR%\stoptime.txt

echo. > %infofile%
echo . ----------------------------- >> %infofile%
echo .         host: %COMPUTERNAME% >> %infofile%
echo .        start: %startdate% %starttime% >> %infofile%
echo .         stop: %stopdate% %stoptime%  >> %infofile%
echo .    run cases: %DIFF_RUNVV% >> %infofile%
echo .make pictures: %DIFF_MAKEPICS% >> %infofile%
echo .        total: %DIFF_TIME% >> %infofile%
echo . ----------------------------- >> %infofile%

if NOT exist %tosummarydir% goto skip_copyfiles
  echo summary   (local): file://%userprofile%/FDS-SMV/Manuals/SMV_Summary/index.html >> %infofile%
  echo summary (windows): https://googledrive.com/host/0B-W-dkXwdHWNUElBbWpYQTBUejQ/index.html >> %infofile%
  echo summary   (linux): https://googledrive.com/host/0B-W-dkXwdHWNN3N2eG92X2taRFk/index.html >> %infofile%
  copy %fromsummarydir%\index*.html %tosummarydir%  1> Nul 2>&1
  copy %fromsummarydir%\images\*.png %tosummarydir%\images 1> Nul 2>&1
  copy %fromsummarydir%\images2\*.png %tosummarydir%\images2 1> Nul 2>&1
:skip_copyfiles
  

cd %CURDIR%

if exist %emailexe% (
  if %havewarnings% == 0 (
    if %haveerrors% == 0 (
      call %email% %mailToSMV% "smokebot build success on %COMPUTERNAME%! %revision%" %infofile%
    ) else (
      echo "start: %startdate% %starttime% " > %infofile%
      echo " stop: %stopdate% %stoptime% " >> %infofile%
      echo. >> %infofile%
      type %errorlog% >> %infofile%
      call %email% %mailToSMV% "smokebot build failure on %COMPUTERNAME%! %revision%" %infofile%
    )
  ) else (
    if %haveerrors% == 0 (
      echo "start: %startdate% %starttime% " > %infofile%
      echo " stop: %stopdate% %stoptime% " >> %infofile%
      echo. >> %infofile%
      type %warninglog% >> %infofile%
      %email% %mailToSMV% "smokebot build success with warnings on %COMPUTERNAME% %revision%" %infofile%
    ) else (
      echo "start: %startdate% %starttime% " > %infofile%
      echo " stop: %stopdate% %stoptime% " >> %infofile%
      echo. >> %infofile%
      type %errorlog% >> %infofile%
      echo. >> %infofile%
      type %warninglog% >> %infofile%
      call %email% %mailToSMV% "smokebot build failure on %COMPUTERNAME%! %revision%" %infofile%
    )
  )
)

echo smokebot_win completed
cd %CURDIR%
pause
exit

:output_abort_message
  echo "***Fatal error: smokebot build failure on %COMPUTERNAME% %revision%"
  if %havemail% == 1 (
    call %email% %mailToSMV% "smokebot build failure on %COMPUTERNAME% %revision%" %errorlog%
  )
exit /b

:: -------------------------------------------------------------
:GET_TIME
:: -------------------------------------------------------------

%gettimeexe% > time.txt
set /p current_time=<time.txt
exit /b 0

:: -------------------------------------------------------------
:GET_DURATION
:: -------------------------------------------------------------
set /a difftime=%3 - %2
set /a diff_h= %difftime%/3600
set /a diff_m= (%difftime% %% 3600 )/60
set /a diff_s= %difftime% %% 60
if %difftime% GEQ 3600 set duration= %diff_h%h %diff_m%m %diff_s%s
if %difftime% LSS 3600 if %difftime% GEQ 60 set duration= %diff_m%m %diff_s%s
if %difftime% LSS 3600 if %difftime% LSS 60 set duration= %diff_s%s
echo %1: %duration% >> %stagestatus%
exit /b 0

:: -------------------------------------------------------------
:is_file_installed
:: -------------------------------------------------------------

  set program=%1
  %program% -help 1>> %OUTDIR%\stage_exist.txt 2>&1
  type %OUTDIR%\stage_exist.txt | find /i /c "not recognized" > %OUTDIR%\stage_count.txt
  set /p nothave=<%OUTDIR%\stage_count.txt
  if %nothave% == 1 (
    echo "***Fatal error: %program% not present"
    echo "***Fatal error: %program% not present" > %errorlog%
    echo "smokebot run aborted"
    call :output_abort_message
    exit /b 1
  )
  exit /b 0

:: -------------------------------------------------------------
  :does_file_exist
:: -------------------------------------------------------------

set file=%1
set outputfile=%2

if NOT exist %file% (
  echo ***fatal error: problem building %file%. Aborting smokebot
  type %outputfile% >> %errorlog%
  call :output_abort_message
  exit /b 1
)
exit /b 0

:: -------------------------------------------------------------
  :find_smokeview_warnings
:: -------------------------------------------------------------

set search_string=%1
set search_file=%2
set stage=%3

findstr /I %search_string% %search_file% | find /V "commands for target" > %OUTDIR%\stage_warning.txt
type %OUTDIR%\stage_warning.txt | find /v /c "kdkwokwdokwd"> %OUTDIR%\stage_nwarning.txt
set /p nwarnings=<%OUTDIR%\stage_nwarning.txt
if %nwarnings% GTR 0 (
  echo %stage% warnings >> %warninglog%
  echo. >> %warninglog%
  type %OUTDIR%\stage_warning.txt >> %warninglog%
  set havewarnings=1
)
exit /b

:: -------------------------------------------------------------
  :find_fds_warnings
:: -------------------------------------------------------------

set search_string=%1
set search_file=%2
set stage=%3

findstr /I %search_string% %search_file% | find /V "mpif.h"  > %OUTDIR%\stage_warning.txt
type %OUTDIR%\stage_warning.txt | find /c ":"> %OUTDIR%\stage_nwarning.txt
set /p nwarnings=<%OUTDIR%\stage_nwarning.txt
if %nwarnings% GTR 0 (
  echo %stage% warnings >> %warninglog%
  echo. >> %warninglog%
  type %OUTDIR%\stage_warning.txt >> %warninglog%
  set havewarnings=1
)
exit /b

:: -------------------------------------------------------------
 :build_guide
:: -------------------------------------------------------------

set guide=%1
set guide_dir=%2

set guideout=%OUTDIR%\stage6_%guide%.txt

cd %guide_dir%

pdflatex -interaction nonstopmode %guide% 1> %guideout% 2>&1
bibtex %guide% 1> %guideout% 2>&1
pdflatex -interaction nonstopmode %guide% 1> %guideout% 2>&1
pdflatex -interaction nonstopmode %guide% 1> %guideout% 2>&1
bibtex %guide% 1>> %guideout% 2>&1

type %guideout% | find "Undefined control" > %OUTDIR%\stage_error.txt
type %guideout% | find "! LaTeX Error:" >> %OUTDIR%\stage_error.txt
type %guideout% | find "Fatal error" >> %OUTDIR%\stage_error.txt
type %guideout% | find "Error:" >> %OUTDIR%\stage_error.txt

type %OUTDIR%\stage_error.txt | find /v /c "JDIJWIDJIQ"> %OUTDIR%\stage_nerrors.txt
set /p nerrors=<%OUTDIR%\stage_nerrors.txt
if %nerrors% GTR 0 (
  echo Errors from Stage 6 - Build %guide% >> %errorlog%
  type %OUTDIR%\stage_error.txt >> %errorlog%
  set haveerrors=1
)

type %guideout% | find "undefined" > %OUTDIR%\stage_warning.txt
type %guideout% | find "multiply"  >> %OUTDIR%\stage_warning.txt

type %OUTDIR%\stage_warning.txt | find /c ":"> %OUTDIR%\nwarnings.txt
set /p nwarnings=<%OUTDIR%\nwarnings.txt
if %nwarnings% GTR 0 (
  echo Warnings from Stage 6 - Build %guide% >> %warninglog%
  type %OUTDIR%\stage_warning.txt >> %warninglog%
  set havewarnings=1
)

copy %guide%.pdf %fromsummarydir%\manuals
copy %guide%.pdf %tosummarydir%\manuals

exit /b

:eof
cd %CURDIR%
pause
exit
