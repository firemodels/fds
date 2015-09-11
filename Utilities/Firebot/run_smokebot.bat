@echo off
:: usage: 
::  run_smokebot -cfastrepo name -fdsrepo name -altemail -email address -nomatlab -noupdate
::  (all command arguments are optional)

set altemail=0
set update=1
set stopscript=0

set cfastrepo=%userprofile%\cfastgitclean
if x%CFASTGIT% == x goto skip_cfastgit
  if EXIST %CFASTGIT% (
    set cfastrepo=%CFASTGIT%
  )
:skip_cfastgit

set fdsrepo=%userprofile%\FDS-SMVgitclean
if exist ..\..\FDS-SMV_root.txt (
  set fdsrepo=..\..
)
if x%FDSGIT% == x goto skip_fdsgit
  if EXIST %FDSGIT% (
    set fdsrepo=%FDSGIT%
  )
:skip_fdsgit

set emailto=
if not x%EMAILGIT% == x (
  set emailto=%EMAILGIT%
)

:: parse command line arguments

call :normalise %cfastrepo% 
set cfastrepo=%temparg%

if %fdsrepo% == none goto skip_fdsrepo
  call :normalise %fdsrepo%
  set fdsrepo=%temparg%
)
:skip_fdsrepo

set stopscript=0
call :getopts %*
if %stopscript% == 1 (
  exit /b
)

:: normalize directory paths

call :normalise %CD% curdir
set curdir=%temparg%

call :normalise %fdsrepo%\Utilities\Firebot
set fdsbotdir=%temparg%

call :normalise %cfastrepo% 
set cfastrepo=%temparg%

if %fdsrepo% == none goto skip_fdsrepo2
  call :normalise %fdsrepo%
  set fdsrepo=%temparg%
:skip_fdsrepo2

set running=%curdir%\bot.running

if exist %running% goto skip_running

:: get latest smokebot

    if %update% == 0 goto no_update
    echo getting latest smokebot
    cd %fdsrepo%
    git fetch origin
    git pull 1> Nul 2>&1
    if not %fdsbotdir% == %curdir% (
      copy %fdsbotdir%\smokebot_win.bat %curdir%
    )
    cd %curdir%
    :no_update

:: run smokebot

  echo 1 > %running%
  call smokebot_win.bat %cfastrepo% %fdsrepo% %update% %altemail% %emailto%
  erase %running%
  goto end_running
:skip_running
  echo smokebot is currently running.
  echo If this is not the case, erase the file %running%
:end_running

goto eof

:getopts
 if (%1)==() exit /b
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
 if /I "%1" EQU "-email" (
   set emailto=%2
   set valid=1
   shift
 )
 if /I "%1" EQU "-altemail" (
   set altemail=1
 )
 if /I "%1" EQU "-noupdate" (
   set valid=1
   set update=0
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
echo -help           - display this message
echo -cfastrepo name - specify the cfast repository
echo       (default: %cfastrepo%) 
echo -fdsrepo name   - specify the FDS-SMV repository
echo       (default: %fdsrepo%) 
echo -altmail        - use an alternate email server
echo -email address  - override "to" email addresses specified in repo 
if "%emailto%" NEQ "" (
echo       (default: %emailto%^)
)
echo -noupdate       - do not update repository
exit /b

:normalise
set temparg=%~f1
exit /b

:eof

