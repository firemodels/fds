@echo off

Rem get FDS-SMV test cases for Windows

set envfile=c:\bin\fds_smv_env.bat
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
cd %svn_root%\Utilities\Scripts\to_google\

set testdir=fds_test_cases_%test_cases_revision%
if exist Test_cases rmdir /s /q Test_cases
svn export https://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Test_cases Test_cases
if exist %testdir%.zip erase %testdir%.zip
wzzip -a -r -P %testdir%.zip Test_cases
if exist %testdir%.exe erase %testdir%.exe
d:\bin\winzip\wzipse32 %testdir%.zip -d "c:\program files\nist\"
erase %testdir%.zip
pause
