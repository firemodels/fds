@echo off

Rem Windows batch file for creating Smokeview User guide figures

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

echo Creating figures for the Smokeview User's guide

call %envfile%

%svn_drive%
cd %svn_root%\Manuals\SMV_User_Guide

pdflatex SMV_User_Guide
bibtex SMV_User_Guide
pdflatex SMV_User_Guide
pdflatex SMV_User_Guide

cd %svn_root%\Manuals\SMV_Verification_Guide

pdflatex SMV_Verification_Guide
bibtex SMV_Verification_Guide
pdflatex SMV_Verification_Guide
pdflatex SMV_Verification_Guide
