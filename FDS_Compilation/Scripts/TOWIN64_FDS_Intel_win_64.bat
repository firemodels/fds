@echo off
Title Building FDS for 32 bit Windows

Rem Batch file used to build a 32 bit version of FDS

set envfile="%homedrive%\%homepath%\fds_smv_env.bat"
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV_5/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%
pause

%svn_drive%
cd %svn_root%\FDS_Compilation\intel_win_64
copy fds5_win_64.exe "c:\program files\fds\fds5\bin\fds5_win_64.exe"
cd %svn_root%\FDS_Compilation\mpi_intel_win_64
copy fds5_mpi_win_64.exe "c:\program files\fds\fds5\bin\fds5_mpi_win_64.exe"
pause
