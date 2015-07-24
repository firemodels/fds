@echo off

set curdir=%CD%
set running=bot.running
set altemail=0
set usematlab=1
set stopscript=0
set fdsrepo=FDS-SMVgitclean
set emailto=

call :getopts %*
if %stopscript% == 1 (
  exit /b
)
  
set gitrepo=%userprofile%\%fdsrepo%
set running=%curdir%\bot.running

if not exist %running% (
  echo 1 > %running%
  call firebot_win.bat %fdsrepo% %altemail% %usematlab% %emailto%
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
 if /I "%1" EQU "-fdsrepo" (
   set fdsrepo=%2
   set valid=1
   shift
 )
 if /I "%1" EQU "-email" (
   set emailto=%2
   set valid=1
   shift
 )
 if /I "%1" EQU "-altemail" (
   set valid=1
   set altemail=1
 )
 if /I "%1" EQU "-nomatlab" (
   set valid=1
   set usematlab=0
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
if not (%1)==() goto getopts
exit /b

:usage  
echo run_smokebot [options]
echo. 
echo -help          - display this message
echo -fdsrepo name  - specify the FDS repo name (default: FDS-SMVgitclean) 
echo -altemail      - use alternate email server
echo -email address - override "to" email addresses specified in repo 
echo -nomatlab      - do not use matlab
exit /b

:eof


