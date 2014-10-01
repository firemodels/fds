@echo off

Rem  Windows batch file to create an achive for a 64 bit Linux smokeview

Rem setup environment variables (defining where repository resides etc) 

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

echo.
echo ---downloading guides
echo.
set manualdir=/home2/smokevis2/firebot/FDS-SMV/Manuals
set firelocaldir="%userprofile%"\FIRE-LOCAL\reports\fds_manuals

pscp %svn_logon%:%manualdir%/FDS_User_Guide/FDS_User_Guide.pdf  %firelocaldir%\.
pscp %svn_logon%:%manualdir%/FDS_Verification_Guide/FDS_Verification_Guide.pdf  %firelocaldir%\.
pscp %svn_logon%:%manualdir%/FDS_Technical_Reference_Guide/FDS_Technical_Reference_Guide.pdf  %firelocaldir%\.
pscp %svn_logon%:%manualdir%/FDS_Validation_Guide/FDS_Validation_Guide.pdf %firelocaldir%\.
pscp %svn_logon%:%manualdir%/FDS_Configuration_Management_Plan/FDS_Configuration_Management_Plan.pdf %firelocaldir%\.

pause
