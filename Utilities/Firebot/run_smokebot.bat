@echo off
set curdir=%CD%
set running=bot.running
set altemail=0
set stopscript=0
set cfastrepo=cfastgitclean
set fdsrepo=FDS-SMVgitclean

call :getopts %*
if %stopscript% == 1 (
  exit /b
)
  
set gitrepo=%userprofile%\%fdsrepo%
set curdir=%CD%
set running=%curdir%\bot.running


if not exist %running% (
  cd %gitrepo%
  git fetch origin
  git pull

  copy Utilities\Firebot\smokebot_win.bat %curdir%
  cd %curdir%

  echo 1 > %running%
  call smokebot_win.bat %cfastrepo% %fdsrepo% %altemail%
  cd %curdir%
  erase %running%
) else (
  echo A bot is already running.
  echo Erase the file %running% if this is not the case
)

goto eof

:getopts
 set valid=0
 set arg=%1
 if /I "%1" EQU "-help" (
   call :usage
   set stopscript=1
   exit /b
 )
 if /I "%1" EQU "-cfastrepo" (
   set cfastrepo=%2
   set valid=1
   shift
 )
 if /I "%1" EQU "-fdsrepo" (
   set fdsrepo=%2
   set valid=1
   shift
 )
 if /I "%1" EQU "-altemail" (
   set altemail=1
 )
 shift
 if %valid% == 0 (
   echo.
   echo ***Error: the input argument %arg% is invalid
   echo.
   echo Usage:
   call :usage
   set stopscript=1
   exit /b
 )
if not (%1)==() goto GETOPTS
exit /b

:usage  
echo run_smokebot [options]
echo. 
echo -help           - display this message
echo -cfastrepo name - specify the cfast repo name (default: cfastgitclean) 
echo -fdsrepo name   - specify the FDS repo name (default: FDS-SMVgitclean) 
echo -email address  - override "to" email addresses specified in repo 
exit /b

:eof


