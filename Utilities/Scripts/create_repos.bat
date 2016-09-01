@echo off
set UN=%username%

set UN=gforney
echo You are about to download the repos:
echo cfast, cor, exp, fds, out, radcal, smv from
echo git@github.com:%UN%
echo Press any key to continue, CTRL c to abort or
echo %0 -h
echo for other options
set /p varname=""

git clone git@github.com:%UN%/cfast.git
git clone git@github.com:%UN%/cor.git
git clone git@github.com:%UN%/exp.git
git clone --recursive git@github.com:%UN%/fds.git
git clone git@github.com:%UN%/out.git
git clone git@github.com:%UN%/radcal.git
git clone git@github.com:%UN%/smv.git

goto eof

:getopts
 if (%1)==() exit /b
 set valid=0
 set arg=%1
 if /I "%1" EQU "-u" (
   set UN=%2
   set valid=1
   shift
 )
 if /I "%1" EQU "-h" (
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
echo Download firemodels repos from github
echo 
echo Options:
echo -u user - specify user [default: %UN%]
echo           specify firemodels to download central repos
echo -h - display this message
exit /b

:eof
