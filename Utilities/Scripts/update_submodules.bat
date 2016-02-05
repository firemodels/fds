@echo off
:: usage: 
::  update_submodules -repo reponame -help
::  (all command arguments are optional)

set curdir=%CD%
set fdsrepo=%userprofile%\FDS-SMVgitclean
if x%FDSGIT% == x goto skip_fdsgit
  if EXIST %FDSGIT% (
    set fdsrepo=%FDSGIT%
  )
:skip_fdsgit

cd ..\..
if EXIST .gitmodules (
  set fdsrepo=%CD%
)
cd %curdir%
  

set stopscript=0
call :getopts %*
if %stopscript% == 1 (
  exit /b
)

cd %fdsrepo%
echo updating submodules in %fdsrepo%
git submodule foreach git remote update 
git submodule foreach git merge origin/master
cd %curdir%

goto eof

:getopts
 if (%1)==() exit /b
 set valid=0
 set arg=%1
 if /I "%1" EQU "-repo" (
   set fdsrepo=%2
   set valid=1
   shift
 )
 if /I "%1" EQU "-help" (
   call :usage
   set stopscript=1
   exit /b
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
echo update_submodules [options]
echo update submodules contained in a git repo.
echo This batch file assumes that the command
echo.
echo git submodule update --init --recursive
echo.
echo was run initially to set up the submodule.
echo. 
echo -help           - display this message
echo -repo name   - specify the git repository
echo       (default: %fdsrepo%) 
exit /b

:normalise
set temparg=%~f1
exit /b

:eof

