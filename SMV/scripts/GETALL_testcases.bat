@echo off

Rem get test cases for Windows and Linux/OSX 

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%
%svn_drive%
cd %svn_root%\Utilities

set LREPOS=FDS-SMV

set testdir=verification_%verification_revision%

Rem get test cases for Linux/OSX systems

plink %svn_logon% %LREPOS%/Utilities/Scripts/get_testcases.csh %verification_revision% %LREPOS%
pscp  %svn_logon%:%LREPOS%/Utilities/uploads/%testdir%.tar.gz %svn_root%\Utilities\uploads\.

Rem get test cases for Windows systems

call %svn_root%\Utilities\Scripts\get_testcases

pause
