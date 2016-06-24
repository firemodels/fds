@echo off

set size=%1

set SCRIPT_DIR=%CD%

cd %CD%\..
set BASEDIR=%CD%

cd %BASEDIR%\..\
set SVNROOT=%CD%
set VDIR=%SVNROOT%\Verification
set INDIR=%SVNROOT%\Verification\Visualization\frames
set OUTDIR=%SVNROOT%\Manuals\SMV_Summary\movies
set WUIINDIR=%SVNROOT%\Verification\WUI\frames
set RUNSMV=call "%SCRIPT_DIR%\runsmv.bat"
set RUNTSMV=call "%SCRIPT_DIR%\runtsmv.bat"
set RUNMSMV=%SCRIPT_DIR%\runmsmv.bat

if "%size%" == "" (
  set SMOKEVIEW=smokeview
  set FDS=fds
) else (
  set SMOKEVIEW=%SVNROOT%\SMV\Build\intel_win_%size%\smokeview_win_%size%.exe -bindir %SVNROOT%\SMV\for_bundle
  set FDS=%SVNROOT%\FDS\Build\intel_win_%size%\fds_win_%size%.exe
)

call :is_file_installed %SMOKEVIEW%|| exit /b 1
call :is_fds_installed %FDS%|| exit /b 1

:: create version string

cd %VDIR%\Visualization
%FDS% version2.fds

cd %VDIR%
call %RUNSMV% Visualization version2

:: -------- plume5c movie -------------------

cd %VDIR%

:: generate movie frames
call %RUNMSMV% Visualization plume5c

cd %INDIR%

:: make movies out of frames generated above

call :frame2movie plume5c_tslice plume5c_tslice

call :frame2movie plume5c_3dsmoke plume5c_3dsmoke

call :frame2movie plume5c_vtslice plume5c_vtslice

call :frame2movie plume5c_iso plume5c_iso

call :frame2movie plume5c_tbound plume5c_tbound

call :frame2movie plume5c_part plume5c_part

:: -------- thouse5 movies -------------------

cd %VDIR%

:: generate movie frames

call %RUNMSMV% Visualization thouse5

cd $INDIR

:: make movies out of frames generated above

call :frame2movie thouse5_tslice thouse5_tslice

call :frame2movie thouse5_smoke3d thouse5_smoke3d

:: -------- BT10m_2x2km_LS movie -------------------

call :make_movie WUI BT10m_2x2km_LS %WUIINDIR%

:: -------- hill_structure movie -------------------

call :make_movie WUI hill_structure %WUIINDIR%

:: -------- levelset1 movie -------------------

call :make_movie WUI levelset1 %WUIINDIR%

:: -------- wind_test1 movie -------------------

call :make_movie WUI wind_test1 %WUIINDIR%

:: -------- tree_test2 movie -------------------

call :make_movie WUI tree_test2 %WUIINDIR%

goto eof

:: -----------------------------------------
  :is_file_installed
:: -----------------------------------------

  set program=%1
  %program% -help 1>> %SCRIPT_DIR%\exist.txt 2>&1
  type %SCRIPT_DIR%\exist.txt | find /i /c "not recognized" > %SCRIPT_DIR%\count.txt
  set /p nothave=<%SCRIPT_DIR%\count.txt
  if %nothave% == 1 (
    echo "***Fatal error: %program% not present"
    echo "Verification aborted"
    exit /b 1
  )
  echo %program% exists
  exit /b 0

:: -----------------------------------------
  :is_fds_installed
:: -----------------------------------------

  set program=%1
  echo "" | %program%  1>> %SCRIPT_DIR%\exist.txt 2>&1
  type %SCRIPT_DIR%\exist.txt | find /i /c "not recognized" > %SCRIPT_DIR%\count.txt
  set /p nothave=<%SCRIPT_DIR%\count.txt
  if %nothave% == 1 (
    echo "***Fatal error: %program% not present"
    echo "Verification aborted"
    exit /b 1
  )
  echo %program% exists
  exit /b 0

:: -----------------------------------------
  :make_movie
:: -----------------------------------------

  set CASEDIR=%1
  set BASEFILE=%2
  set FRAMEDIR=%3

  cd %VDIR%

:: generate movie frames

  call %RUNMSMV% %CASEDIR% %BASEFILE%

  cd %FRAMEDIR%

:: make movies out of frames generated above

  call :frame2movie %BASEFILE%_movie %BASEFILE%_movie

  exit /b 0

:: -----------------------------------------
  :frame2movie
:: -----------------------------------------

  set BASEFILE=%1
  set BASEMOVIE=%2
  ffmpeg -y -f image2 -i "%BASEFILE%_%%04d.png" %OUTDIR%\%BASEMOVIE%.m1v
  exit /b 0

:eof
