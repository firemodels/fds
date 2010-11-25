@echo off

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

echo test prerequisites used by the FDS/Smokeview web scripts
echo press any key to continue
pause >NUL

%svn_drive%
call %envfile%

echo.
echo Press any key to show configuration settings
pause >NUL
call %svn_root%\SMV\scripts\SHOW_setup.bat

echo.
echo Testing putty installation: running plink -V
pause
plink -V

echo.
echo Testing Compiler installation
pause
IF DEFINED IFORT_COMPILER12 echo compiler located at %IFORT_COMPILER12%
IF NOT DEFINED IFORT_COMPILER12 echo IFORT_COMPILER12 variable not defined

IF DEFINED IFORT_COMPILER11 echo compiler located at %IFORT_COMPILER11%
IF NOT DEFINED IFORT_COMPILER11 echo IFORT_COMPILER11 variable not defined

echo.
echo Testing Winzip installation
echo First run wzzip
wzzip
pause
echo Second run winzipse
winzipse

echo.
echo Testing complete
pause

