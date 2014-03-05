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
echo ---downloading images
echo.
erase %svn_root%\Manuals\SMV_Summary\images\*.png
pscp %svn_logon%:FDS-SMV/Manuals/SMV_Summary/images/*.png  %svn_root%\Manuals\SMV_Summary\images\.

echo.
echo ---downloading animations
echo.
erase %svn_root%\Manuals\SMV_Summary\movies\*.m1v
erase %svn_root%\Manuals\SMV_Summary\movies\*.png
pscp %svn_logon%:FDS-SMV/Manuals/SMV_Summary/movies/*.m1v  %svn_root%\Manuals\SMV_Summary\movies\.
pscp %svn_logon%:FDS-SMV/Manuals/SMV_Summary/movies/*.png  %svn_root%\Manuals\SMV_Summary\movies\.

erase %svn_root%\Manuals\SMV_Summary\movies2\*.m1v
erase %svn_root%\Manuals\SMV_Summary\movies2\*.png
pscp %svn_logon%:FDS-SMV/Manuals/SMV_Summary/movies2/*.m1v  %svn_root%\Manuals\SMV_Summary\movies2\.
pscp %svn_logon%:FDS-SMV/Manuals/SMV_Summary/movies2/*.png  %svn_root%\Manuals\SMV_Summary\movies2\.

%svn_drive%
cd %svn_root%\smv\scripts
set version=%smv_version%

set scriptdir=FDS-SMV/SMV/scripts
set UG=Manuals\SMV_User_Guide
set VG=Manuals\SMV_Verification_Guide
set TG=Manuals\SMV_Technical_Reference_Guide

echo.
echo ---downloading guides
echo.
pscp %svn_logon%:FDS-SMV/Manuals/SMV_User_Guide/SMV_User_Guide.pdf  %svn_root%\%UG%\.
pscp %svn_logon%:FDS-SMV/Manuals/SMV_Verification_Guide/SMV_Verification_Guide.pdf  %svn_root%\%VG%\.
pscp %svn_logon%:FDS-SMV/Manuals/SMV_Technical_Reference_Guide/SMV_Technical_Reference_Guide.pdf %svn_root%\%TG%\.

pause
