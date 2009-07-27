@echo off

Rem get FDS-SMV test cases for Windows

set envfile="%homedrive%\%homepath%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV_5/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

Rem --------- should not need to edit below ----------

call %envfile%
%svn_drive%
cd %svn_root%\Utilities\to_google\
set excludefile=%svn_root%\Utilities\Scripts\examples_win.exclude
set testdir=verification_%verification_revision%

if exist Examples rmdir /s /q Examples
svn export https://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Verification Examples
cd Examples
if exist ..\%testdir%.zip erase ..\%testdir%.zip
wzzip -a -r -P ..\%testdir%.zip *
if exist ..\%testdir%.exe erase ..\%testdir%.exe
Rem c:\bin\winzip\wzipse32 ..\%testdir%.zip -d "c:\program files\nist\Examples"
c:\bin\winzip\wzipse32 ..\%testdir%.zip -d "c:\program files\fds5\Examples"
erase ..\%testdir%.zip
pause
