@echo off
set prog=%1

:: batch file to build smokeview utility programs on Windows, Linux or OSX platforms

:: setup environment variables (defining where repository resides etc) 

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use smv/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%
echo.
echo  Building %prog% for 64 bit %platform%
Title Building %prog% for 64 bit %platform%

%svn_drive%


cd %svn_root%\smv\Build\%prog%\intel_win_64
call make_%prog%

plink %linux_logon% %linux_svn_root%/smv/scripts/run_command.sh smv/Build/%prog%/intel_linux_64 make_%prog%.sh

plink %osx_logon% %linux_svn_root%/smv/scripts/run_command.sh smv/Build/%prog%/intel_osx_64 make_%prog%.sh
