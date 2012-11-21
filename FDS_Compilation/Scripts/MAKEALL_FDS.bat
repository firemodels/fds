@echo off
Title Building Parallel FDS for 32 bit Windows

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
set winscriptdir=%svn_root%\FDS_Compilation\Scripts

call %winscriptdir%\CLEAN_FDS_Intel_win_32

call %winscriptdir%\CLEAN_FDS_Intel_linux_32

call %winscriptdir%\CLEAN_FDS_Intel_linux_64

call %winscriptdir%\CLEAN_FDS_Intel_OSX_32

call %winscriptdir%\CLEAN_FDS_Intel_OSX_64

echo "Building FDS for all platforms (except for 64 bit windows)

call %winscriptdir%\MAKE_FDS_Intel_win_32
call %winscriptdir%\MAKE_FDS_Intel_win_mpi32

call %winscriptdir%\MAKE_FDS_Intel_linux_32
call %winscriptdir%\MAKE_FDS_Intel_linux_mpi32

call %winscriptdir%\MAKE_FDS_Intel_linux_64
call %winscriptdir%\MAKE_FDS_Intel_linux_mpi64

call %winscriptdir%\MAKE_FDS_Intel_OSX_32
call %winscriptdir%\MAKE_FDS_Intel_OSX_mpi32

call %winscriptdir%\MAKE_FDS_Intel_OSX_64
