@echo off
Title Cleaning FDS for 32 bit Windows 

Rem Batch file used to clean 32 and 64 bit FDS build directories

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
echo cleaning All FDS for all platforms
call %svn_root%\FDS_Compilation\Scripts\CLEAN_FDS_Intel_win_32.bat
call %svn_root%\FDS_Compilation\Scripts\CLEAN_FDS_Intel_win_64.bat
call %svn_root%\FDS_Compilation\Scripts\CLEAN_FDS_Intel_linux_32.bat
call %svn_root%\FDS_Compilation\Scripts\CLEAN_FDS_Intel_linux_64.bat
call %svn_root%\FDS_Compilation\Scripts\CLEAN_FDS_Intel_OSX_32.bat
call %svn_root%\FDS_Compilation\Scripts\CLEAN_FDS_Intel_OSX_64.bat

pause