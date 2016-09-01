@echo off

echo You are about to setup remotes for the repos:
echo cfast, cor, exp, fds, out, radcal, smv at 
echo git@github.com:firemodels
echo Press any key to continue or CTRL c to abort.
echo You should only run this command if you have cloned these
echo repos from repos that have been forked from firemodels
set /p varname=""

CURDIR=%CD%

cd %CD%\cor
git remote add firemodels git@github.com:firemodels/cfast.git
git remote update

git remote add firemodels git@github.com:firemodels/cor.git
git remote update

cd %CD%\exp
git remote add firemodels git@github.com:firemodels/exp.git
git remote update

cd %CD%\fds
git remote add firemodels git@github.com:firemodels/fds.git
git remote update

cd %CD%\out
git remote add firemodels git@github.com:firemodels/out.git
git remote update

cd %CD%\radcal
git remote add firemodels git@github.com:firemodels/radcal.git
git remote update

cd %CD%\smv
git remote add firemodels git@github.com:firemodels/smv.git
git remote update

cd %CURDIR%
goto eof

:getopts
 if (%1)==() exit /b
 set valid=0
 set arg=%1
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
echo Setup remote connections with each firemodels repo
echo.
echo Options:
echo -h - display this message
exit /b
