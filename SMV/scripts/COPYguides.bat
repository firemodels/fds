@echo off
setlocal EnableDelayedExpansion
set whichguides=%1

::  batch file to copy guides

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

%svn_drive%


if "%whichguides%" == "from_smv_local" (
  set fromdir=%svn_root%\Manuals
  set todir="%userprofile%"\FDS_Guides

  Title Copying smokeview guides from local repo to %userprofile%\FDS_Guides
  echo copying guides
  echo from directory: !fromdir!
  echo to directory: !todir!

  call :COPY !fromdir!\SMV_User_Guide\SMV_User_Guide.pdf                                !todir!\.
  call :COPY !fromdir!\SMV_Verification_Guide\SMV_Verification_Guide.pdf                !todir!\.
  call :COPY !fromdir!\SMV_Technical_Reference_Guide\SMV_Technical_Reference_Guide.pdf  !todir!\.
  goto eof
)
if "%whichguides%" == "from_smv_linux" (
  set fromdir=%smokebotrepo%/Manuals
  set todir="%userprofile%"\FDS_Guides

  Title Copying smokeview guides from linux smokebot
  echo copying guides
  echo from directory: !fromdir!
  echo to directory: !todir!

  pscp %linux_logon%:!fromdir!/SMV_User_Guide/SMV_User_Guide.pdf                                !todir!\.
  pscp %linux_logon%:!fromdir!/SMV_Verification_Guide/SMV_Verification_Guide.pdf                !todir!\.
  pscp %linux_logon%:!fromdir!/SMV_Technical_Reference_Guide/SMV_Technical_Reference_Guide.pdf  !todir!\.
  goto eof
)
if "%whichguides%" == "from_fds_linux" (
  set fromdir=%firebotrepo%/Manuals
  set todir="%userprofile%"\FDS_Guides

  Title Copying FDS guides from linux firebot
  echo copying guides
  echo from directory: !fromdir!
  echo to directory: !todir!

  pscp %linux_logon%:!fromdir!/FDS_Config_Management_Plan/FDS_Config_Management_Plan.pdf        !todir!\.
  pscp %linux_logon%:!fromdir!/FDS_Technical_Reference_Guide/FDS_Technical_Reference_Guide.pdf  !todir!\.
  pscp %linux_logon%:!fromdir!/FDS_User_Guide/FDS_User_Guide.pdf                                !todir!\.
  pscp %linux_logon%:!fromdir!/FDS_Verification_Guide/FDS_Verification_Guide.pdf                !todir!\.
  pscp %linux_logon%:!fromdir!/FDS_Validation_Guide/FDS_Validation_Guide.pdf                    !todir!\.
  goto eof
)
if "%whichguides%" == "to_linux" (
  set fromdir="%userprofile%"\FDS_Guides
  set todir=FDS_Guides

  Title uploading guides to Linux:%linux_hostname%
  echo uploading guides to %linux_hostname%
  echo from directory: !fromdir!
  echo to directory: !todir!

  pscp !fromdir!\*.pdf %linux_logon%:!todir!/.
  goto eof
)
if "%whichguides%" == "to_osx" (
  set fromdir="%userprofile%"\FDS_Guides
  set todir=FDS_Guides

  Title uploading guides to OSX:%osx_hostname%
  echo uploading guides to %osx_hostname%
  echo from directory: !fromdir!
  echo to directory: !todir!

  pscp !fromdir!\*.pdf %osx_logon%:!todir!/.
  goto eof
)
goto eof

:COPY
set label=%~n1%~x1
set infile=%1
set infiletime=%~t1
set outfile=%2
IF EXIST %infile% (
   echo copying %label% %infiletime%
   copy %infile% %outfile% >Nul
) ELSE (
   echo.
   echo *** warning: %infile% does not exist
   echo.
   pause
)
exit /b


:eof
echo.
echo copy complete
pause
