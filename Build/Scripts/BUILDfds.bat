@echo off
set platform=%1
set buildtype=%2

:: batch file to build FDS on Windows, Linux or OSX platforms

:: setup environment variables (defining where repository resides etc) 

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
echo.
echo  Building %buildtype% FDS for %platform%
Title Building %buildtype% FDS for %platform%

%svn_drive%

set type=
if "%buildtype%" == "debug" (
  set mpi=
  set type=_db
)
if "%buildtype%" == "dev" (
  set mpi=
  set type=_dv
)
if "%buildtype%" == "release" (
  set mpi=mpi_
  set type=
)

if "%platform%" == "windows" (
  cd %svn_root%\FDS_Compilation\%mpi%intel_win_64%type%
  if "%buildtype%" == "release" (
    erase *.obj *.mod *.exe
  )
  call make_fds
  goto eof
)
if "%platform%" == "linux" (
  if "%buildtype%" == "release" (
    plink %linux_logon% %linux_svn_root%/SMV/scripts/run_command.sh FDS_Compilation/%mpi%intel_linux_64%type% clean_fds.sh
  )
  plink %linux_logon% %linux_svn_root%/SMV/scripts/run_command.sh FDS_Compilation/%mpi%intel_linux_64%type% make_fds.sh
  pause
  goto eof
)
if "%platform%" == "osx" (
  if "%buildtype%" == "release" (
    plink %osx_logon% %linux_svn_root%/SMV/scripts/run_command.sh FDS_Compilation/%mpi%intel_osx_64%type% clean_fds.sh
  )
  plink %osx_logon% %linux_svn_root%/SMV/scripts/run_command.sh FDS_Compilation/%mpi%intel_osx_64%type% make_fds.sh
  pause
  goto eof
)

:eof
