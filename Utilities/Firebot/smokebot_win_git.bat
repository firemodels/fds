@echo   off

::  set number of OpenMP threads

set OMP_NUM_THREADS=1

:: -------------------------------------------------------------
::                         set repository names
:: -------------------------------------------------------------

set fdsbasename=FDS-SMVclean
set svnroot=%userprofile%\%fdsbasename%
if NOT exist %svnroot% (
  cd %userprofile%
  echo %fdsbasename% repository does not exist - creating.
  svn co http://fds-smv.googlecode.com/svn/trunk/FDS/trunk %fdsbasename%
  echo %fdsbasename% repository creation complete.
)

set cfastbasename=cfastclean
set cfastroot=%userprofile%\%cfastbasename%
if NOT exist %cfastroot% (
  cd %userprofile%
  echo %cfastbasename% repository does not exist - creating.
  svn co http://cfast.googlecode.com/svn/trunk/cfast/trunk %cfastbasename%
  echo %cfastbasename% repository creation complete.
)

:: -------------------------------------------------------------
::                         setup environment
:: -------------------------------------------------------------

set CURDIR=%CD%

if not exist output mkdir output
if not exist history mkdir history
if not exist timings mkdir timings

set OUTDIR=%CURDIR%\output
set HISTORYDIR=%CURDIR%\history
set TIMINGSDIR=%CURDIR%\timings
set timefile=%OUTDIR%\time.txt

erase %OUTDIR%\*.txt 1> Nul 2>&1

set email=%svnroot%\SMV\scripts\email.bat

set debug=1
set release=0
set errorlog=%OUTDIR%\stage_errors.txt
set warninglog=%OUTDIR%\stage_warnings.txt
set errorwarninglog=%OUTDIR%\stage_errorswarnings.txt
set infofile=%OUTDIR%\stage_info.txt
set revisionfilestring=%OUTDIR%\revision.txt
set revisionfilenum=%OUTDIR%\revision_num.txt
set stagestatus=%OUTDIR%\stage_status.log

set fromsummarydir=%svnroot%\Manuals\SMV_Summary
set tosummarydir="%SMOKEBOT_SUMMARY_DIR%"

set haveerrors=0
set havewarnings=0
set haveCC=1

set emailexe=%userprofile%\bin\mailsend.exe
set gettimeexe=%svnroot%\Utilities\get_time\intel_win_64\get_time.exe
set runbatchexe=%svnroot%\SMV\source\runbatch\intel_win_64\runbatch.exe

date /t > %OUTDIR%\starttime.txt
set /p startdate=<%OUTDIR%\starttime.txt
time /t > %OUTDIR%\starttime.txt
set /p starttime=<%OUTDIR%\starttime.txt

call "%svnroot%\Utilities\Scripts\setup_intel_compilers.bat" 1> Nul 2>&1
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

call :GET_TIME TIME_beg
call :GET_TIME PRELIM_beg

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

call :is_file_installed grep|| exit /b 1
echo             found grep

call :is_file_installed sed|| exit /b 1
echo             found sed

call :is_file_installed cut|| exit /b 1
echo             found cut

call :is_file_installed svn|| exit /b 1
echo             found svn

echo. 1> %OUTDIR%\stage0.txt 2>&1

:: revert cfast repository

cd %cfastroot%
if "%cfastbasename%" == "cfastclean" (
   echo             reverting %cfastbasename% repository
   cd %cfastroot%
   call :svn_revert 1>> %OUTDIR%\stage0.txt 2>&1
)

:: update cfast repository

echo             updating %cfastbasename% repository
cd %cfastroot%
svn update  1>> %OUTDIR%\stage0.txt 2>&1

:: revert FDS/Smokeview repository

cd %svnroot%
if "%fdsbasename%" == "FDS-SMVclean" (
   echo             reverting %fdsbasename% repository
   cd %svnroot%
   call :svn_revert 1>> %OUTDIR%\stage0.txt 2>&1
)

:: update FDS/Smokeview repository

echo             updating %fdsbasename% repository
svn update 1>> %OUTDIR%\stage0.txt 2>&1

svn info | grep Revision > %revisionfilestring%
set /p revisionstring=<%revisionfilestring%

svn info | grep Revision | cut -d " " -f 2 > %revisionfilenum%
set /p revisionnum=<%revisionfilenum%

set errorlogpc=%HISTORYDIR%\errors_%revisionnum%.txt
set warninglogpc=%HISTORYDIR%\warnings_%revisionnum%.txt

set timingslogfile=%TIMINGSDIR%\timings_%revisionnum%.txt

:: build cfast

echo             building cfast
cd %cfastroot%\CFAST\intel_win_64
erase *.obj *.mod *.exe 1>> %OUTDIR%\stage0.txt 2>&1
make VPATH="../Source:../Include" INCLUDE="../Include" -f ..\makefile intel_win_64 1>> %OUTDIR%\stage0.txt 2>&1
call :does_file_exist cfast7_win_64.exe %OUTDIR%\stage0.txt|| exit /b 1

call :GET_DURATION PRELIM %PRELIM_beg%

:: -------------------------------------------------------------
::                           stage 1
:: -------------------------------------------------------------

call :GET_TIME BUILDFDS_beg

echo Stage 1 - Building FDS

echo             parallel debug

cd %svnroot%\FDS_Compilation\mpi_intel_win_64_db
erase *.obj *.mod *.exe 1> %OUTDIR%\stage1b.txt 2>&1
make VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_64_db 1>> %OUTDIR%\stage1b.txt 2>&1

call :does_file_exist fds_mpi_win_64_db.exe %OUTDIR%\stage1b.txt|| exit /b 1
call :find_fds_warnings "warning" %OUTDIR%\stage1b.txt "Stage 1b"

echo             parallel release

cd %svnroot%\FDS_Compilation\mpi_intel_win_64
erase *.obj *.mod *.exe 1> %OUTDIR%\stage1d.txt 2>&1
make VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_64  1>> %OUTDIR%\stage1d.txt 2>&1

call :does_file_exist fds_mpi_win_64.exe %OUTDIR%\stage1d.txt|| exit /b 1
call :find_fds_warnings "warning" %OUTDIR%\stage1d.txt "Stage 1d"

call :GET_DURATION BUILDFDS %BUILDFDS_beg%

:: -------------------------------------------------------------
::                           stage 2
:: -------------------------------------------------------------

call :GET_TIME BUILDSMVUTIL_beg

echo Stage 2 - Building Smokeview

echo             libs

cd %svnroot%\SMV\Build\LIBS\lib_win_intel_64
call makelibs2 1>> %OUTDIR%\stage2a.txt 2>&1

echo             debug

cd %svnroot%\SMV\Build\intel_win_64
erase *.obj *.mod *.exe smokeview_win_64_db.exe 1> %OUTDIR%\stage2a.txt 2>&1
make -f ..\Makefile intel_win_64_db 1>> %OUTDIR%\stage2a.txt 2>&1

call :does_file_exist smokeview_win_64_db.exe %OUTDIR%\stage2a.txt|| exit /b 1
call :find_smokeview_warnings "warning" %OUTDIR%\stage2a.txt "Stage 2a"

echo             release

cd %svnroot%\SMV\Build\intel_win_64
erase *.obj *.mod smokeview_win_64.exe 1> %OUTDIR%\stage2b.txt 2>&1
make -f ..\Makefile intel_win_64 1>> %OUTDIR%\stage2b.txt 2>&1

call :does_file_exist smokeview_win_64.exe %OUTDIR%\stage2b.txt|| aexit /b 1
call :find_smokeview_warnings "warning" %OUTDIR%\stage2b.txt "Stage 2b"

:: -------------------------------------------------------------
::                           stage 3
:: -------------------------------------------------------------

echo Stage 3 - Building FDS/Smokeview utilities

echo             fds2ascii
cd %svnroot%\Utilities\fds2ascii\intel_win_64
erase *.obj *.mod *.exe 1> %OUTDIR%\stage3c.txt 2>&1
ifort -o fds2ascii_win_64.exe /nologo ..\..\Data_processing\fds2ascii.f90  1>> %OUTDIR%\stage3.txt 2>&1
call :does_file_exist fds2ascii_win_64.exe %OUTDIR%\stage3.txt|| exit /b 1

if %haveCC% == 1 (
  echo             background
  cd %svnroot%\Utilities\background\intel_win_32
  erase *.obj *.mod *.exe 1>> %OUTDIR%\stage3.txt 2>&1
  make -f ..\Makefile intel_win_32 1>> %OUTDIR%\stage3.txt 2>&1
  call :does_file_exist background.exe %OUTDIR%\stage3.txt

  echo             smokediff
  cd %svnroot%\Utilities\smokediff\intel_win_64
  erase *.obj *.mod *.exe 1>> %OUTDIR%\stage3.txt 2>&1
  make -f ..\Makefile intel_win_64 1>> %OUTDIR%\stage3.txt 2>&1
  call :does_file_exist smokediff_win_64.exe %OUTDIR%\stage3.txt

  echo             smokezip
  cd %svnroot%\Utilities\smokezip\intel_win_64
  erase *.obj *.mod *.exe 1>> %OUTDIR%\stage3.txt 2>&1
  make -f ..\Makefile intel_win_64 1>> %OUTDIR%\stage3.txt 2>&1
  call :does_file_exist smokezip_win_64.exe %OUTDIR%\stage3.txt|| exit /b 1

  echo             wind2fds
  cd %svnroot%\Utilities\wind2fds\intel_win_64
  erase *.obj *.mod *.exe 1>> %OUTDIR%\stage3.txt 2>&1
  make -f ..\Makefile intel_win_64 1>> %OUTDIR%\stage3.txt 2>&1
  call :does_file_exist wind2fds_win_64.exe %OUTDIR%\stage3.txt|| exit /b 1
) else (
  call :is_file_installed background|| exit /b 1
  echo             background not built, using installed version
  call :is_file_installed smokediff|| exit /b 1
  echo             smokediff not built, using installed version
  call :is_file_installed smokezip|| exit /b 1
  echo             smokezip not built, using installed version
  call :is_file_installed wind2fds|| exit /b 1
  echo             wind2fds not built, using installed version
)

call :GET_DURATION PRELIM %PRELIM_beg%

:: -------------------------------------------------------------
::                           stage 4
:: -------------------------------------------------------------

call :GET_TIME RUNVV_beg

echo Stage 4 - Running verification cases
echo             debug mode

:: run the cases

cd %svnroot%\Verification\scripts
call Run_SMV_cases %debug% 1> %OUTDIR%\stage4a.txt 2>&1

:: check the cases

cd %svnroot%\Verification\scripts
echo. > %OUTDIR%\stage_error.txt
call Check_SMV_cases 

:: report errors

call :report_errors Stage 4a, "Debug FDS case errors"|| exit /b 1

echo             release mode

:: run the cases

cd %svnroot%\Verification\scripts
call Run_SMV_cases %release% 1> %OUTDIR%\stage4b.txt 2>&1

:: check the cases

cd %svnroot%\Verification\scripts
echo. > %OUTDIR%\stage_error.txt
call Check_SMV_cases 

:: report errors

call :report_errors Stage 4b, "Release FDS case errors"|| exit /b 1

call :GET_DURATION RUNVV %RUNVV_beg%

:: -------------------------------------------------------------
::                           stage 5
:: -------------------------------------------------------------

call :GET_TIME MAKEPICS_beg

echo Stage 5 - Making Smokeview pictures

cd %svnroot%\Verification\scripts
call MAKE_SMV_pictures 64 1> %OUTDIR%\stage5.txt 2>&1

call :find_smokeview_warnings "error" %OUTDIR%\stage5.txt "Stage 5"

call :GET_DURATION MAKEPICS %MAKEPICS_beg%

:: -------------------------------------------------------------
::                           stage 6
:: -------------------------------------------------------------

call :GET_TIME MAKEGUIDES_beg

echo Stage 6 - Building Smokeview guides

echo             Technical Reference
call :build_guide SMV_Technical_Reference_Guide %svnroot%\Manuals\SMV_Technical_Reference_Guide 1>> %OUTDIR%\stage6.txt 2>&1

echo             Verification
call :build_guide SMV_Verification_Guide %svnroot%\Manuals\SMV_Verification_Guide 1>> %OUTDIR%\stage6.txt 2>&1

echo             User
call :build_guide SMV_User_Guide %svnroot%\Manuals\SMV_User_Guide 1>> %OUTDIR%\stage6.txt 2>&1

echo             Geom Notes
call :build_guide geom_notes %svnroot%\Manuals\FDS_User_Guide 1>> %OUTDIR%\stage6.txt 2>&1

call :GET_DURATION MAKEGUIDES %MAKEGUIDES_beg%
call :GET_DURATION TOTALTIME %TIME_beg%

:: -------------------------------------------------------------
::                           wrap up
:: -------------------------------------------------------------

date /t > %OUTDIR%\stoptime.txt
set /p stopdate=<%OUTDIR%\stoptime.txt
time /t > %OUTDIR%\stoptime.txt
set /p stoptime=<%OUTDIR%\stoptime.txt

echo. > %infofile%
echo . -----------------------------         >> %infofile%
echo .         host: %COMPUTERNAME%          >> %infofile%
echo .        start: %startdate% %starttime% >> %infofile%
echo .         stop: %stopdate% %stoptime%   >> %infofile%
echo .        setup: %DIFF_PRELIM%           >> %infofile%
echo .    run cases: %DIFF_RUNVV%            >> %infofile%
echo .make pictures: %DIFF_MAKEPICS%         >> %infofile%
echo .  make guides: %DIFF_MAKEGUIDES%       >> %infofile%
echo .        total: %DIFF_TOTALTIME%        >> %infofile%
echo . -----------------------------         >> %infofile%

copy %infofile% %timingslogfile%

if NOT exist %tosummarydir% goto skip_copyfiles
  echo summary   (local): file://%userprofile%/FDS-SMV/Manuals/SMV_Summary/index.html >> %infofile%
  echo summary (windows): https://googledrive.com/host/0B-W-dkXwdHWNUElBbWpYQTBUejQ/index.html >> %infofile%
  echo summary   (linux): https://googledrive.com/host/0B-W-dkXwdHWNN3N2eG92X2taRFk/index.html >> %infofile%
  copy %fromsummarydir%\index*.html %tosummarydir%  1> Nul 2>&1
  copy %fromsummarydir%\images\*.png %tosummarydir%\images 1> Nul 2>&1
  copy %fromsummarydir%\images2\*.png %tosummarydir%\images2 1> Nul 2>&1
:skip_copyfiles
  

cd %CURDIR%

sed "s/$/\r/" < %warninglog% > %warninglogpc%
sed "s/$/\r/" < %errorlog% > %errorlogpc%

if exist %emailexe% (
  if %havewarnings% == 0 (
    if %haveerrors% == 0 (
      call %email% %mailToSMV% "smokebot success on %COMPUTERNAME%! %revisionstring%" %infofile%
    ) else (
      echo "start: %startdate% %starttime% " > %infofile%
      echo " stop: %stopdate% %stoptime% " >> %infofile%
      echo. >> %infofile%
      type %errorlogpc% >> %infofile%
      call %email% %mailToSMV% "smokebot failure on %COMPUTERNAME%! %revisionstring%" %infofile%
    )
  ) else (
    if %haveerrors% == 0 (
      echo "start: %startdate% %starttime% " > %infofile%
      echo " stop: %stopdate% %stoptime% " >> %infofile%
      echo. >> %infofile%
      type %warninglogpc% >> %infofile%
      %email% %mailToSMV% "smokebot success with warnings on %COMPUTERNAME% %revisionstring%" %infofile%
    ) else (
      echo "start: %startdate% %starttime% " > %infofile%
      echo " stop: %stopdate% %stoptime% " >> %infofile%
      echo. >> %infofile%
      type %errorlogpc% >> %infofile%
      echo. >> %infofile%
      type %warninglogpc% >> %infofile%
      call %email% %mailToSMV% "smokebot failure on %COMPUTERNAME%! %revisionstring%" %infofile%
    )
  )
)

echo smokebot_win completed
goto :eof

:output_abort_message
  echo "***Fatal error: smokebot failure on %COMPUTERNAME% %revisionstring%"
  if exist %emailexe% (
    call %email% %mailToSMV% "smokebot failure on %COMPUTERNAME% %revisionstring%" %errorlog%
  )
exit /b

:: -------------------------------------------------------------
:report_errors
:: -------------------------------------------------------------
set stage_label=%1
grep -v " " %OUTDIR%\stage_error.txt | wc -l > %OUTDIR%\stage_nerror.txt
set /p nerrors=<%OUTDIR%\stage_nerror.txt
if %nerrors% GTR 0 (
   echo %stage_label% >> %errorlog%
   echo. >> %errorlog%
   type %OUTDIR%\stage_error.txt >> %errorlog%
   set haveerrors=1
   set haveerrors_now=1
   call :output_abort_message
   exit /b 1
)
exit /b 0

:: -------------------------------------------------------------
:GET_DURATION
:: -------------------------------------------------------------

:: compute difftime=time2 - time1

set label=%1
set time1=%2

set difftime=DIFF_%label%
call :GET_TIME time2

set /a diff=%time2% - %time1%
set /a diff_h= %diff%/3600
set /a diff_m= (%diff% %% 3600 )/60
set /a diff_s= %diff% %% 60
if %diff% GEQ 3600 set duration= %diff_h%h %diff_m%m %diff_s%s
if %diff% LSS 3600 if %diff% GEQ 60 set duration= %diff_m%m %diff_s%s
if %diff% LSS 3600 if %diff% LSS 60 set duration= %diff_s%s
echo %label%: %duration% >> %stagestatus%
set %difftime%=%duration%
exit /b 0

:: -------------------------------------------------------------
:GET_TIME
:: -------------------------------------------------------------

set arg1=%1

%gettimeexe% > %timefile%
set /p %arg1%=<%timefile%
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

grep -v "commands for target" %search_file% > %OUTDIR%\stage_warning0.txt
grep -i -A 5 -B 5 %search_string% %OUTDIR%\stage_warning0.txt > %OUTDIR%\stage_warning.txt
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

grep -v "mpif.h" %search_file% > %OUTDIR%\stage_warning0.txt
grep -i -A 5 -B 5 %search_string% %OUTDIR%\stage_warning0.txt  > %OUTDIR%\stage_warning.txt
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
:svn_revert
:: -------------------------------------------------------------
svn cleanup .
svn revert -R .
For /f "tokens=1,2" %%A in ('svn status --no-ignore') Do (
     If [%%A]==[?] ( Call :UniDelete %%B
     ) Else If [%%A]==[I] Call :UniDelete %%B
   )
exit /b

:: -------------------------------------------------------------
:UniDelete delete file/dir
:: -------------------------------------------------------------
if "%1"=="%~nx0" exit /b
IF EXIST "%1\*" ( 
    RD /S /Q "%1"
) Else (
    If EXIST "%1" DEL /S /F /Q "%1"
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
