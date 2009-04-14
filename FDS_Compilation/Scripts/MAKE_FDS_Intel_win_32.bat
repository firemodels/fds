@echo off
Title Building FDS for 32 bit Windows

Rem Batch file used to build a 32 bit version of FDS

set envfile=%homedrive%\%homepath%\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV_5/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

Rem location of batch files used to set up Intel compilation environment

set intelbin=c:\bin

Rem call %intelbin%\iclvars ia32
call %intelbin%\clvars x86
call %intelbin%\ifortvars ia32

call %envfile%

%svn_drive%
cd %svn_root%\FDS_Compilation\intel_win_32
set out=intel_win_32.out
echo. | tee -a %out%
date /t | tee -a  %out%
time /t | tee -a  %out%
make VPATH="../../FDS_Source" -f ..\makefile intel_win_32 | tee -a %out%

pause
