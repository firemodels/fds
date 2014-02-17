@echo off

set size=%1

::ffmpeg -f image2 -i "thouse5_%%04d.png" mpg_video.mpg
::ffmpeg -f image2 -i "thouse5_%%04d.png" -intra -coder ac -an mpgh_video.mpg

set SCRIPT_DIR=%CD%

cd %CD%\..
set BASEDIR=%CD%

cd %BASEDIR%\..\
set SVNROOT=%CD%
set VDIR=%SVNROOT%\Verification
set INDIR=%SVNROOT%\Verification\Visualization\frames
set OUTDIR=%SVNROOT%\Manuals\SMV_Summary\movies
set WUIINDIR=%SVNROOT%\Verification\WUI\frames

set RUNMSMV=%SCRIPT_DIR%\runmsmv.bat

if "%size%" == "" (
  set SMOKEVIEW=smokeview
) else (
  set SMOKEVIEW=%SVNROOT%\SMV\Build\intel_win_%size%\smokeview_win_%size%.exe -bindir %SVNROOT%\SMV\for_bundle
)

call :is_file_installed %SMOKEVIEW%|| exit /b 1

set RUNSMV=call "%SCRIPT_DIR%\runsmv.bat"



:: -------- BT10m_2x2km_LS movie -------------------

:: call :make_movie WUI BT10m_2x2km_LS %WUIINDIR%

:: -------- thouse5 movie -------------------

call :make_movie Visualization thouse5 %INDIR%

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
  :make_movie
:: -----------------------------------------

  set CASEDIR=%1
  set BASEFILE=%2
  set FRAMEDIR=%3

  cd %VDIR%

:: generate movie frames
  echo generating movie frames
  call %RUNMSMV% %CASEDIR% %BASEFILE%

  cd %FRAMEDIR%

:: make movies out of frames generated above
  echo making %BASEFILE% movie
  ffmpeg -f image2 -i "%BASEFILE%_%%04d.png" %OUTDIR%\%BASEFILE%_movie.m1v

  exit /b 0

:eof
