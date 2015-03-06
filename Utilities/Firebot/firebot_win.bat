@echo off

:: $Date$ 
:: $Revision$
:: $Author$

set DEBUGOPT=%1

:: -------------------------------------------------------------
::                         set environment
:: -------------------------------------------------------------

:: set number of OpenMP threads

set OMP_NUM_THREADS=1

:: -------------------------------------------------------------
::                         set repository names
:: -------------------------------------------------------------

set fdsbasename=FDS-SMVclean
set svnroot=%userprofile%\%fdsbasename%
if NOT exist %svnroot% (
  echo ***Fatal error: The svn repository %fdsbasename% does not exist.
  echo Aborting firebot
  exit /b 1
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

erase %OUTDIR%\*.txt %OUTDIR%\*.log 1> Nul 2>&1

set email=%svnroot%\SMV\scripts\email.bat

set release=0
set debug=1
set errorlog=%OUTDIR%\firebot_errors.txt
set timefile=%OUTDIR%\time.txt
set datefile=%OUTDIR%\date.txt
set waitfile=%OUTDIR%\wait.txt
set warninglog=%OUTDIR%\firebot_warnings.txt
set errorwarninglog=%OUTDIR%\firebot_errorswarnings.txt
set infofile=%OUTDIR%\firebot_info.txt
set revisionfilestring=%OUTDIR%\revision_string.txt
set revisionfilenum=%OUTDIR%\revision_num.txt
set stagestatus=%OUTDIR%\firebot_status.log
set counta=%OUTDIR%\firebot_count0a.txt
set countb=%OUTDIR%\firebot_count0b.txt
set scratchfile=%OUTDIR%\firebot_scratch.txt

set fromsummarydir=%svnroot%\Manuals\SMV_Summary

set haveerrors=0
set havewarnings=0
set have_icc=1

set emailexe=%userprofile%\bin\mailsend.exe
set gettimeexe=%svnroot%\Utilities\get_time\intel_win_64\get_time.exe
set runbatchexe=%svnroot%\SMV\source\runbatch\intel_win_64\runbatch.exe

call :get_datetime startdate starttime

call "%svnroot%\Utilities\Scripts\setup_intel_compilers.bat" 1> Nul 2>&1
call %svnroot%\Utilities\Firebot\firebot_email_list.bat
if "%DEBUGOPT%" == "debug" (
   set mailToList=%mailToFDSDebug%
) else (
   set mailToList=%mailToFDS%
)

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

call :get_time TIME_beg
call :get_time PRELIM_beg

ifort 1> %scratchfile% 2>&1
type %scratchfile% | find /i /c "not recognized" > %counta%
set /p nothave_ifort=<%counta%
set have_ifort=1
if %nothave_ifort% == 1 (
  echo "***Fatal error: Fortran compiler not present"
  echo "***Fatal error: Fortran compiler not present" > %errorlog%
  echo "firebot run aborted"
  call :output_abort_message
  exit /b 1
)
echo             found Fortran

icl 1> %scratchfile% 2>&1
type %scratchfile% | find /i /c "not recognized" > %countb%
set /p nothave_icc=<%countb%
if %nothave_icc% == 1 (
  set have_icc=0
  echo "***Warning: C/C++ compiler not found - using installed Smokeview to generate images"
) else (
  echo             found C/C++
)

if NOT exist %emailexe% (
  echo ***Warning: email client not found.   
  echo             firebot messages will only be sent to the console.
) else (
  echo             found mailsend
)

call :is_file_installed pdflatex|| exit /b 1
echo             found pdflatex

call :is_file_installed grep|| exit /b 1
echo             found grep

call :is_file_installed sed|| exit /b 1
echo             found sed

call :is_file_installed svn|| exit /b 1
echo             found svn

echo. 1>> %OUTDIR%\stage0.txt 2>&1

:: revert FDS/Smokeview repository

if "%fdsbasename%" == "FDS-SMVclean" (
   echo             reverting %fdsbasename% repository
   cd %svnroot%
   call :svn_revert 1> Nul 2>&1
)

:: update FDS/Smokeview repository

echo             updating %fdsbasename% repository
cd %svnroot%
svn update 1>> %OUTDIR%\stage0.txt 2>&1

svn info | grep Revision > %revisionfilestring%
set /p revisionstring=<%revisionfilestring%

svn info | grep Revision | cut -d " " -f 2 > %revisionfilenum%
set /p revisionnum=<%revisionfilenum%

set errorlogpc=%HISTORYDIR%\errors_%revisionnum%.txt
set warninglogpc=%HISTORYDIR%\warnings_%revisionnum%.txt

set timingslogfile=%TIMINGSDIR%\timings_%revisionnum%.txt

call :get_time PRELIM_end
call :get_duration PRELIM DIFF_PRELIM %PRELIM_end% %PRELIM_beg%

:: -------------------------------------------------------------
::                           stage 1
:: -------------------------------------------------------------

call :get_time BUILDFDS_beg

echo Stage 1 - Building FDS

echo             parallel debug

cd %svnroot%\FDS_Compilation\mpi_intel_win_64_db
erase *.obj *.mod *.exe *.pdb 1> Nul 2>&1
make VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_64_db 1> %OUTDIR%\makefdsd.log 2>&1
call :does_file_exist fds_mpi_win_64_db.exe %OUTDIR%\makefdsd.log|| exit /b 1
call :find_warnings "warning" %OUTDIR%\makefdsd.log "Stage 1b, FDS parallel debug compilation"

echo             parallel release

cd %svnroot%\FDS_Compilation\mpi_intel_win_64
erase *.obj *.mod *.exe *.pdb 1> Nul 2>&1
make VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_64  1> %OUTDIR%\makefdsr.log 2>&1
call :does_file_exist fds_mpi_win_64.exe %OUTDIR%\makefdsr.log|| exit /b 1
call :find_warnings "warning" %OUTDIR%\makefdsr.log "Stage 1d, FDS parallel release compilation"

call :get_time BUILDFDS_end
call :get_duration BUILDFDS DIFF_BUILDFDS %BUILDFDS_end% %BUILDFDS_beg%

:: -------------------------------------------------------------
::                           stage 2
:: -------------------------------------------------------------

call :get_time BUILDSMVUTIL_beg

echo Stage 2 - Building Smokeview

echo             libs

cd %svnroot%\SMV\Build\LIBS\lib_win_intel_64
call makelibs2 1>> %OUTDIR%\stage2a.txt 2>&1

echo             debug

cd %svnroot%\SMV\Build\intel_win_64
erase *.obj *.mod *.exe smokeview_win_64_db.exe 1> Nul 2>&1
make -f ..\Makefile intel_win_64_db 1> %OUTDIR%\makesmvd.log 2>&1
call :does_file_exist smokeview_win_64_db.exe %OUTDIR%\makesmvd.log|| exit /b 1
call :find_warnings "warning" %OUTDIR%\makesmvd.log "Stage 2a, Smokeview debug compilation"

echo             release

cd %svnroot%\SMV\Build\intel_win_64
erase *.obj *.mod smokeview_win_64.exe 1> Nul 2>&1
make -f ..\Makefile intel_win_64 1> %OUTDIR%\makesmvr.log 2>&1

call :does_file_exist smokeview_win_64.exe %OUTDIR%\makesmvr.log|| aexit /b 1
call :find_warnings "warning" %OUTDIR%\makesmvr.log "Stage 2b, Smokeview release compilation"

:: -------------------------------------------------------------
::                           stage 3
:: -------------------------------------------------------------

echo Stage 3 - Building Utilities

echo             fds2ascii
cd %svnroot%\Utilities\fds2ascii\intel_win_64
erase *.obj *.mod *.exe 1> Nul 2>&1
ifort -o fds2ascii_win_64.exe /nologo ..\..\Data_processing\fds2ascii.f90  1> %OUTDIR%\makefds2ascii.log 2>&1
call :does_file_exist fds2ascii_win_64.exe %OUTDIR%\makefds2ascii.log|| exit /b 1
call :find_warnings "warning" %OUTDIR%\makefds2ascii.log "Stage 3, Building FDS/Smokeview utilities"

if %have_icc% == 1 (
  echo             background
  cd %svnroot%\Utilities\background\intel_win_32
  erase *.obj *.mod *.exe 1> Nul 2>&1
  make -f ..\Makefile intel_win_32 1> %OUTDIR%\makebackground.log 2>&1
  call :does_file_exist background.exe %OUTDIR%\makebackground.log
  call :find_warnings "warning" %OUTDIR%\makebackground.log "Stage 3, Building FDS/Smokeview utilities"
) else (
  call :is_file_installed background|| exit /b 1
  echo             background not built, using installed version
)

call :get_time BUILDSMVUTIL_end=%current_time% 
call :get_duration BUILDSMVUTIL DIFF_BUILDSMVUTIL %BUILDSMVUTIL_end% %BUILDSMVUTIL_beg%

:: -------------------------------------------------------------
::                           stage 4
:: -------------------------------------------------------------

call :get_time RUNVV_beg

echo Stage 4 - Running verification cases
echo             debug mode

:: run cases

cd %svnroot%\Verification\
call Run_FDS_cases %debug% 1> %OUTDIR%\stage4a.txt 2>&1

:: check cases

set haveerrors_now=0
echo. > %OUTDIR%\stage_error.txt
cd %svnroot%\Verification\
call Check_FDS_cases 

:: report errors

call :report_errors Stage 4a, "Debug FDS case errors"|| exit /b 1

echo             release mode

:: run cases

cd %svnroot%\Verification\
call Run_FDS_cases %release% 1> %OUTDIR%\stage4b.txt 2>&1

:: check cases

set haveerrors_now=0
echo. > %OUTDIR%\stage_error.txt
cd %svnroot%\Verification\
call Check_FDS_cases

:: report errors

call :report_errors Stage 4b, "Release FDS case errors"|| exit /b 1

call :get_time RUNVV_end
call :get_duration RUNVV DIFF_RUNVV %RUNVV_end% %RUNVV_beg%

:: -------------------------------------------------------------
::                           stage 5
:: -------------------------------------------------------------

if exist %emailexe% (
  call :get_datetime current_ddate current_ttime

  echo "making pictures on %COMPUTERNAME% %revisionstring%" > %infofile%
  echo "  start time: %startdate% %starttime%" >> %infofile%
  echo "current time: %current_ddate% %current_ttime%" >> %infofile%
  call %email% %mailToFDSDebug% "making pictures on %COMPUTERNAME% %revisionstring%" %infofile%
)
call :get_time MAKEPICS_beg

echo Stage 5 - Making pictures
echo             FDS verification cases

cd %svnroot%\Verification\
call MAKE_FDS_pictures 64 1> %OUTDIR%\stage5.txt 2>&1

call :get_time MAKEPICS_end
call :get_duration MAKEPICS DIFF_MAKEPICS %MAKEPICS_end% %MAKEPICS_beg%

:: -------------------------------------------------------------
::                           stage 6
:: -------------------------------------------------------------

call :get_time MAKEGUIDES_beg

:: echo Stage 6 - Building guides

:: don't build FDS guides until a "matlab" stage is added
::echo             FDS Technical Reference
::call :build_guide FDS_Technical_Reference_Guide %svnroot%\Manuals\FDS_Technical_Reference_Guide 1>> %OUTDIR%\stage6.txt 2>&1

::echo             FDS User
::call :build_guide FDS_User_Guide %svnroot%\Manuals\FDS_User_Guide 1>> %OUTDIR%\stage6.txt 2>&1

::echo             FDS Verification
::call :build_guide FDS_Verification_Guide %svnroot%\Manuals\FDS_Verification_Guide 1>> %OUTDIR%\stage6.txt 2>&1

::echo             FDS Validation
::call :build_guide FDS_Validation_Guide %svnroot%\Manuals\FDS_Validation_Guide 1>> %OUTDIR%\stage6.txt 2>&1

call :get_time MAKEGUIDES_end
call :get_duration MAKEGUIDES DIFF_MAKEGUIDES %MAKEGUIDES_end% %MAKEGUIDES_beg%

call :get_time TIME_end
call :get_duration TOTALTIME DIFF_TIME %TIME_end% %TIME_beg%

:: -------------------------------------------------------------
::                           wrap up
:: -------------------------------------------------------------

call :get_datetime stopdate stoptime

echo. > %infofile%
echo . ----------------------------- >> %infofile%
echo .         host: %COMPUTERNAME% >> %infofile%
echo .        start: %startdate% %starttime% >> %infofile%
echo .         stop: %stopdate% %stoptime%  >> %infofile%
echo .    run cases: %DIFF_RUNVV% >> %infofile%
echo .make pictures: %DIFF_MAKEPICS% >> %infofile%
echo .        total: %DIFF_TIME% >> %infofile%
echo . ----------------------------- >> %infofile%

copy %infofile% %timingslogfile%

cd %CURDIR%

sed "s/$/\r/" < %warninglog% > %warninglogpc%
sed "s/$/\r/" < %errorlog% > %errorlogpc%

if exist %emailexe% (
  if %havewarnings% == 0 (
    if %haveerrors% == 0 (
      call %email% %mailToList% "firebot success on %COMPUTERNAME%! %revisionstring%" %infofile%
    ) else (
      echo "start: %startdate% %starttime% " > %infofile%
      echo " stop: %stopdate% %stoptime% " >> %infofile%
      echo. >> %infofile%
      type %errorlogpc% >> %infofile%
      call %email% %mailToList% "firebot failure on %COMPUTERNAME%! %revisionstring%" %infofile%
    )
  ) else (
    if %haveerrors% == 0 (
      echo "start: %startdate% %starttime% " > %infofile%
      echo " stop: %stopdate% %stoptime% " >> %infofile%
      echo. >> %infofile%
      type %warninglogpc% >> %infofile%
      call %email% %mailToList% "firebot success with warnings on %COMPUTERNAME% %revisionstring%" %infofile%
    ) else (
      echo "start: %startdate% %starttime% " > %infofile%
      echo " stop: %stopdate% %stoptime% " >> %infofile%
      echo. >> %infofile%
      type %errorlogpc% >> %infofile%
      echo. >> %infofile%
      type %warninglogpc% >> %infofile%
      call %email% %mailToList% "firebot failure on %COMPUTERNAME%! %revisionstring%" %infofile%
    )
  )
)

echo firebot_win completed
goto :eof

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
:output_abort_message
:: -------------------------------------------------------------
   for %%a in (%errorlogpc%) do (
      set filename=%%~nxa
   )    

  echo ***Fatal error: firebot failure on %COMPUTERNAME% %revisionstring%, error log: %filename%
  sed "s/$/\r/" < %errorlog% > %errorlogpc%
  if exist %emailexe% (
    echo sending email to %mailToList%
    call %email% %mailToList% "firebot failure on %COMPUTERNAME% %revisionstring%" %errorlogpc%
  )
exit /b

:: -------------------------------------------------------------
:get_datetime
:: -------------------------------------------------------------

set arg1=%1
set arg2=%2

date /t > %datefile%
set /p %arg1%=<%datefile%

time /t > %timefile%
set /p %arg2%=<%timefile%

exit /b 0

:: -------------------------------------------------------------
:get_time
:: -------------------------------------------------------------

set arg1=%1

%gettimeexe% > %timefile%
set /p %arg1%=<%timefile%
exit /b 0

:: -------------------------------------------------------------
:get_duration
:: -------------------------------------------------------------

:: compute difftime=time2 - time1

set label=%1
set difftime=%2
set time2=%3
set time1=%4

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
:is_file_installed
:: -------------------------------------------------------------

  set program=%1
  %program% -help 1>> %OUTDIR%\stage_exist.txt 2>&1
  type %OUTDIR%\stage_exist.txt | find /i /c "not recognized" > %OUTDIR%\stage_count.txt
  set /p nothave=<%OUTDIR%\stage_count.txt
  if %nothave% == 1 (
    echo "***Fatal error: %program% not present"
    echo "***Fatal error: %program% not present" > %errorlog%
    echo "firebot run aborted"
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
  echo ***Fatal error: problem building %file%. Aborting firebot
  echo ***Fatal error: problem building %file%. Aborting firebot >> %errorlog%
  type %outputfile% >> %errorlog%
  call :output_abort_message
  exit /b 1
)
exit /b 0

:: -------------------------------------------------------------
:find_warnings
:: -------------------------------------------------------------

set search_string=%1
set search_file=%2
set stage=%3

grep -v "commands for target" %search_file% > %OUTDIR%\stage_warning0.txt
grep -v "mpif.h" %OUTDIR%\stage_warning0.txt > %OUTDIR%\stage_warning1.txt
grep -i -A 5 -B 5 %search_string% %OUTDIR%\stage_warning1.txt > %OUTDIR%\stage_warning.txt
type %OUTDIR%\stage_warning.txt | find /v /c "  "> %OUTDIR%\stage_nwarning.txt
set /p nwarnings=<%OUTDIR%\stage_nwarning.txt
if %nwarnings% GTR 0 (
  echo %stage% warnings >> %warninglog%
  echo. >> %warninglog%
  type %OUTDIR%\stage_warning.txt >> %warninglog%
  set havewarnings=1
)
exit /b

:: -------------------------------------------------------------
:find_errors
:: -------------------------------------------------------------

set search_string=%1
set search_file=%2
set stage=%3

grep -v "commands for target" %search_file% > %OUTDIR%\stage_error0.txt
grep -v "mpif.h" %OUTDIR%\stage_error0.txt > %OUTDIR%\stage_error1.txt
grep -i -A 5 -B 5 %search_string% %OUTDIR%\stage_error1.txt > %OUTDIR%\stage_error.txt
type %OUTDIR%\stage_error.txt | find /v /c "  "> %OUTDIR%\stage_nerror.txt
set /p nerrors=<%OUTDIR%\stage_nerror.txt
if %nerrors% GTR 0 (
  echo %stage% errors >> %errorlog%
  echo. >> %errorlog%
  type %OUTDIR%\stage_error.txt >> %errorlog%
  set haveerrors=1
  set haveerrors_now=1
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

type %OUTDIR%\stage_error.txt | find /v /c "  "> %OUTDIR%\stage_nerrors.txt
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

exit /b

:eof
cd %CURDIR%
