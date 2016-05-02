@echo off

echo Testing that software and environment variables used by 
echo the FDS/Smokeview web scripts are setup.  
echo.
echo Press any key to begin test.
pause >NUL

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo The file fds_smv_env.bat was not found in your
echo home directory: %userprofile%
echo.
echo 1.  Copy it from SMV\for_bundle in the repository to %userprofile%
echo 2.  Then edit it to match your configuration.
echo 3.  Then repeat this test.
echo. Press any key to continue.
pause>NUL
goto:eof

:endif_envexist

echo fds_smv_env.bat was found in your home directory.
echo Press any key to view settings.
pause >NUL

%svn_drive%
call %envfile%

call %svn_root%\SMV\scripts\SHOW_setup.bat

echo.
echo Press any key to test putty installation
pause>NUL
plink -V

echo.
echo Press any key to test Compiler installation
pause>NUL
IF DEFINED IFORT_COMPILER16 echo IFORT_COMPILER16 variable is defined.
IF DEFINED IFORT_COMPILER16 echo Fortran compiler located at %IFORT_COMPILER16%
IF NOT DEFINED IFORT_COMPILER16 echo IFORT_COMPILER16 variable not defined

echo.
echo Winzip test
echo First press any key to test command line winzip (wzzip)
pause>NUL
wzzip

echo.
echo Next press any key to test the self extractor
pause>NUL
winzipse

echo.
echo Testing complete
pause

