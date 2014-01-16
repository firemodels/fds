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

echo copying images
erase %svn_root%\Manuals\SMV_Summary\images\*.png
copy %svn_root%\Manuals\SMV_User_Guide\SCRIPT_FIGURES\*.png          %svn_root%\Manuals\SMV_Summary\images\.
copy %svn_root%\Manuals\SMV_Verification_Guide\SCRIPT_FIGURES\*.png  %svn_root%\Manuals\SMV_Summary\images\.

pause
