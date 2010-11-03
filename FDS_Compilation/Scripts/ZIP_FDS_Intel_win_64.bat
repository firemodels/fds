@echo off

Rem Batch file used to create a self-extracting archive containing FDS

set envfile="%homedrive%\%homepath%"\fds_smv_env.bat
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
cd %svn_root%\FDS_Compilation

set fdsroot=fds_%fds_version%_win64
set togoogle=%svn_root%\FDS_Compilation\to_google\%fdsroot%

IF EXIST %togoogle% rmdir /S /Q %togoogle%
mkdir %togoogle%

copy Intel_Win_64\fds5_win_64.exe %togoogle%\fds5_win64.exe
copy Mpi_Intel_Win_64\fds5_win_mpi_64.exe %togoogle%\fds5_mpi_win64.exe

echo.
echo winzipping distribution directory
cd %togoogle%
if exist ..\%fdsroot%.zip erase ..\%fdsroot%.zip
wzzip -a -r -P ..\%fdsroot%.zip *
cd ..
echo.
echo creating self-extracting archive
if exist %fdsroot%.exe erase %fdsroot%.exe
wzipse32 %fdsroot%.zip -d "C:\Program Files (x86)\nist\FDS"

echo %fdsroot%.exe located in %cd%
pause
