@echo off

Rem Batch file used to create a self-extracting archive containing FDS

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
cd %svn_root%\Utilities\Makefile

set fdsroot=fds_%fds_version%_%fds_revision%_win32
set togoogle=%svn_root%\Utilities\to_google\%fdsroot%
mkdir %togoogle%
copy Intel_Win_32\fds_win_32.exe %togoogle%\fds6.exe
copy Mpi_Intel_Win_32\fds_win_mpi_32.exe %togoogle%\fds6_mpi.exe

echo.
echo winzipping distribution directory
cd %togoogle%
if exist ..\%fdsroot%.zip erase ..\%fdsroot%.zip
wzzip -a -r -P ..\%fdsroot%.zip *
cd ..
echo.
echo creating self-extracting archive
if exist %fdsroot%.exe erase %fdsroot%.exe
wzipse32 %fdsroot%.zip -d "C:\Program Files\nist\FDS"

echo %fdsroot%.exe located in %cd%
pause
