@echo off

:: $Date$ 
:: $Revision$
:: $Author$

set rundebug=%1
if "%rundebug%" == "1" (
  set DEBUG=_db
) else (
  set DEBUG=
)

set BASEDIR="%CD%"
cd ..
set SVNROOT="%CD%"
cd %BASEDIR%
set TIME_FILE="%BASEDIR%\fds_case_times.txt"

::*** uncomment following two lines to use OpenMP

:: set OMP_NUM_THREADS=1

:: default FDS location

set BACKGROUNDEXE=%SVNROOT%\Utilities\background\intel_win_32\background.exe
set FDSBASE=fds_mpi_win_64%DEBUG%.exe
set FDSEXE=%SVNROOT%\FDS_Compilation\mpi_intel_win_64%DEBUG%\%FDSBASE%

call :is_file_installed %BACKGROUNDEXE%|| exit /b 1
call :is_file_installed %FDSEXE%|| exit /b 1

set FDS=%BACKGROUNDEXE% -u 60 -m 70 -d 5 %FDSEXE%
set QFDS=call %SVNROOT%\Utilities\Scripts\runfds.bat

echo.
echo Creating FDS case list from FDS_Cases.sh
..\Utilities\Data_processing\sh2bat FDS_Cases.sh FDS_Cases.bat

echo.
echo Running FDS cases
echo.

echo "FDS test cases begin" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

call FDS_Cases.bat

goto eof

:: -----------------------------------------
:is_file_installed
:: -----------------------------------------
  set program=%1
  %program% -help 1> %BASEDIR%\exist.txt 2>&1
  type %BASEDIR%\exist.txt | find /i /c "not recognized" > %BASEDIR%\count.txt
  set /p nothave=<%BASEDIR%\count.txt
  if %nothave% GTR 0 (
    echo "***Fatal error: %program% not present"
    echo "Verification suite aborted"
    exit /b 1
  )
  exit /b 0


:eof
echo "FDS test cases end" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

:eof2


