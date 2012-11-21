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

set win_makefile=%svn_root%\Utilities\Makefile
set verify=%svn_root%\Utilities\Makefile\Verify_Build
set linux_makefile=%linux_svn_root%/Utilities/Makefile

echo.
echo erasing previous summary files
if exist %verify%\intel_win_32.out erase %verify%\intel_win_32.out
if exist %verify%\intel_win_64.out erase %verify%\intel_win_64.out
if exist %verify%\mpi_intel_win_32.out erase %verify%\mpi_intel_win_32.out
if exist %verify%\mpi_intel_win_64.out erase %verify%\mpi_intel_win_64.out

if exist     %verify%\intel_linux_32.out erase     %verify%\intel_linux_32.out
if exist     %verify%\intel_linux_64.out erase     %verify%\intel_linux_64.out
if exist %verify%\mpi_intel_linux_32.out erase %verify%\mpi_intel_linux_32.out
if exist %verify%\mpi_intel_linux_64.out erase %verify%\mpi_intel_linux_64.out

if exist     %verify%\intel_osx_32.out erase     %verify%\intel_osx_32.out
if exist     %verify%\intel_osx_64.out erase     %verify%\intel_osx_64.out
if exist %verify%\mpi_intel_osx_32.out erase %verify%\mpi_intel_osx_32.out
Rem if exist %verify%\mpi_intel_osx_64.out erase %verify%\mpi_intel_osx_64.out

echo.
echo copying 32 bit Windows summary files
copy /y %win_makefile%\Intel_Win_32\intel_win_32.out %verify%\.
copy /y %win_makefile%\Mpi_Intel_Win_32\mpi_intel_win_32.out %verify%\.

echo.
echo downloading 64 bit Windows summary files
pscp %svn_logon%:%linux_makefile%/Intel_Win_64/intel_win_64.out %verify%\.
pscp %svn_logon%:%linux_makefile%/Mpi_Intel_Win_64/mpi_intel_win_64.out %verify%\.

echo.
echo downloading Linux summary files
pscp %svn_logon%:%linux_makefile%/Intel_Linux_32/intel_linux_32.out %verify%\.
pscp %svn_logon%:%linux_makefile%/Intel_Linux_64/intel_linux_64.out %verify%\.
pscp %svn_logon%:%linux_makefile%/Mpi_Intel_Linux_32/mpi_intel_linux_32.out %verify%\.
pscp %svn_logon%:%linux_makefile%/Mpi_Intel_Linux_64/mpi_intel_linux_64.out %verify%\.

echo.
echo downloading OSX summary files
pscp %svn_logon%:%linux_makefile%/Intel_OSX_32/intel_osx_32.out %verify%\.
pscp %svn_logon%:%linux_makefile%/Intel_OSX_64/intel_osx_64.out %verify%\.
pscp %svn_logon%:%linux_makefile%/Mpi_Intel_OSX_32/mpi_intel_osx_32.out %verify%\.
Rem pscp %svn_logon%:%linux_makefile%/Mpi_Intel_OSX_64/mpi_intel_osx_64.out %verify%\.

echo.
echo assembling all build summaries into one document
cd %verify%
pdflatex fds_verify_build>Nul
pdflatex fds_verify_build>Nul

echo.
echo opening summary document
start acrord32 %verify%\fds_verify_build.pdf


pause