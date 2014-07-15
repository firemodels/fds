@echo off
Title Building Parallel FDS for 64 bit Windows 
	       
Rem Batch file used to build a 32 bit version of FDS

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

cd %svn_root%\FDS_Compilation\impi_intel_win_64

erase *.obj 
erase *.mod
.\make_fds_from_webpage

pause
