@echo off

Rem  Windows batch file to copy smokediffset SVNROOT=~/FDS-SMV

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

echo.
echo --------------------------------------------
echo Building 64 bit Windows versions of wind2fds
echo --------------------------------------------

cd %svn_root%\Utilities\wind2fds\intel_win_64
call "%IFORT_COMPILER14%\bin\compilervars" intel64
erase *.obj *.mod
make -f ..\Makefile intel_win_64

echo.
echo --------------------------------------------
echo Building 32 bit Windows versions of wind2fds
echo --------------------------------------------

cd %svn_root%\Utilities\wind2fds\intel_win_32
call "%IFORT_COMPILER14%\bin\compilervars" ia32
erase *.obj *.mod
make -f ..\Makefile intel_win_32
pause
