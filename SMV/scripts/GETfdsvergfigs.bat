@echo off

set envfile=%userprofile%\fds_smv_env.bat
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

set CURDIR=%CD%

Title Download FDS verification guide images


cd %svn_root%\Manuals\FDS_Verification_Guide\SCRIPT_FIGURES
pscp %linux_logon%:%firebotrepo%/Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/* .

Title Download FDS verification guide scatterplot images

cd %svn_root%\Manuals\FDS_Verification_Guide\SCRIPT_FIGURES\Scatterplots
pscp %linux_logon%:%firebotrepo%/Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/Scatterplots/* .

cd %CURDIR%
pause
