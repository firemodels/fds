@echo off

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
set fromdir=%firebotrepo%/Manuals
set todir="%userprofile%"\FDS_Guides

pscp %linux_logon%:%fromdir%/FDS_User_Guide/FDS_User_Guide.pdf                                %todir%\.
pscp %linux_logon%:%fromdir%/FDS_Verification_Guide/FDS_Verification_Guide.pdf                %todir%\.
pscp %linux_logon%:%fromdir%/FDS_Technical_Reference_Guide/FDS_Technical_Reference_Guide.pdf  %todir%\.
pscp %linux_logon%:%fromdir%/FDS_Validation_Guide/FDS_Validation_Guide.pdf                    %todir%\.
pscp %linux_logon%:%fromdir%/FDS_Config_Management_Plan/FDS_Config_Management_Plan.pdf        %todir%\.

pause
