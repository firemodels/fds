@echo off

Rem get test cases for Windows and Linux/OSX 

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

call %envfile%
cd %svn_root%/Utilities/Scripts

set LREPOS=FDS-SMV

set testdir=FDS_Test_cases_%test_cases_revision%

Rem get test cases for Linux/OSX systems

plink %svn_logon% %LREPOS%/Utilities/Scripts/get_testcases.csh %test_cases_revision%
pscp %svn_logon%:%LREPOS%/Utilities/Scripts/to_google/%testdir%.tar.gz to_google\.


Rem get test cases for Windows systems

call %svn_root%\Utilities\Scripts\get_testcases

pause
