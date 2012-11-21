@echo off
Title Setup environment for 64 bit Intel compiles
	       
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

Rem location of batch files used to set up Intel compilation environment

set intelbin=c:\bin

IF EXIST "C:\Program Files (x86)\Intel\Compiler\11.0" GOTO setup_64
IF EXIST "C:\Program Files\Intel\Compiler\11.0" GOTO setup_64_on_32
echo 64 bit build environment not available
pause
goto exit

:setup_64
echo setting up  64 bit build environment on 64 bit platform
call %intelbin%\iclvars intel64
call %intelbin%\ifortvars intel64
goto exit

:setup_64_on_32
echo setting up  64 bit build environment on 32 bit platform
call %intelbin%\iclvars ia32_intel64
call %intelbin%\ifortvars ia32_intel64
:exit





