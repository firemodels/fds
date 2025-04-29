@echo off
set CURDIR=%CD%

if "x%FDS_SAVE_PATH%" == "x" goto set_fds_path
:: restore original path

echo restoring original path
set PATH=%FDS_SAVE_PATH%
echo.
echo PATH=%PATH%
set FDS_SAVE_PATH=
echo.
echo original path restored

goto eof

:set_fds_path

:: set path containing fds, smokeview, intel mpi and windows entries

SET I_MPI_ROOT=%~dp0\mpi
SET FDS_PATH=%~dp0
SET SMV_PATH=%~dp0\..\..\SMV6

cd %I_MPI_ROOT%
SET I_MPI_ROOT=%CD%

cd %FDS_PATH%
SET FDS_PATH=%CD%

cd %SMV_PATH%
SET SMV_PATH=%CD%

SET FDS_SAVE_PATH=%PATH%
echo setting path to contain only fds, smokeview, intel mpi and windows entries
echo type fdspath again to restore original path
SET PATH=%SMV_PATH%;%FDS_PATH%;%I_MPI_ROOT%;%WINDIR%;%WINDIR%\system32
echo.
echo PATH=%PATH%

:eof

cd %CURDIR%
