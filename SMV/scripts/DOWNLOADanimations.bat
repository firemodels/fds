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

echo downloading animations
erase %svn_root%\Manuals\SMV_Summary\movies\*.m1v
erase %svn_root%\Manuals\SMV_Summary\movies\*.png
pscp %svn_logon%:FDS-SMV/Manuals/SMV_Summary/movies/*.m1v  %svn_root%\Manuals\SMV_Summary\movies\.
pscp %svn_logon%:FDS-SMV/Manuals/SMV_Summary/movies/*.png  %svn_root%\Manuals\SMV_Summary\movies\.

erase %svn_root%\Manuals\SMV_Summary\movies2\*.m1v
erase %svn_root%\Manuals\SMV_Summary\movies2\*.png
pscp %svn_logon%:FDS-SMV/Manuals/SMV_Summary/movies2/*.m1v  %svn_root%\Manuals\SMV_Summary\movies2\.
pscp %svn_logon%:FDS-SMV/Manuals/SMV_Summary/movies2/*.png  %svn_root%\Manuals\SMV_Summary\movies2\.


pause
