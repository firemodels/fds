@echo off

:: setup environment variables (defining where repository resides etc) 

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use smv/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...

pause>NUL
goto:eof

:endif_envexist

call %envfile%

%svn_drive%

set pubdir=%firebotrepo%/fds/Manuals
pscp %linux_hostname%:%pubdir%/FDS_Config_Management_Plan/FDS_Config_Management_Plan.pdf       %userprofile%\FDS_Guides
pscp %linux_hostname%:%pubdir%/FDS_Technical_Reference_Guide/FDS_Technical_Reference_Guide.pdf %userprofile%\FDS_Guides
pscp %linux_hostname%:%pubdir%/FDS_User_Reference_Guide/FDS_User_Reference_Guide.pdf           %userprofile%\FDS_Guides
pscp %linux_hostname%:%pubdir%/FDS_Validation_Guide/FDS_Validation_Guide.pdf                   %userprofile%\FDS_Guides
pscp %linux_hostname%:%pubdir%/FDS_Verification_Guide/FDS_Verification_Guide.pdf               %userprofile%\FDS_Guides
