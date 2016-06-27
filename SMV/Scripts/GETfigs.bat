@echo off
setlocal EnableDelayedExpansion
set whichguides=%1

::  batch to copy smokview/smokebot or fdsfirebot figures to local repo

::  setup environment variables (defining where repository resides etc) 

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV/Scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%
echo.

%svn_drive%

if "%whichguides%" == "smvug" (
  Title Download smokeview user guide images

  cd %svn_root%\Manuals\SMV_User_Guide\SCRIPT_FIGURES
  pscp %linux_logon%:%smokebotrepo%/Manuals/SMV_User_Guide/SCRIPT_FIGURES/* .
  goto eof
)
if "%whichguides%" == "smvvg" (
  Title Download smokeview verification guide images

  cd %svn_root%\Manuals\SMV_Verification_Guide\SCRIPT_FIGURES
  pscp %linux_logon%:%smokebotrepo%/Manuals/SMV_Verification_Guide/SCRIPT_FIGURES/* .
  goto eof
)
if "%whichguides%" == "fdsug" (
  Title Download FDS user guide images

  cd %svn_root%\Manuals\FDS_User_Guide\SCRIPT_FIGURES
  pscp %linux_logon%:%firebotrepo%/Manuals/FDS_User_Guide/SCRIPT_FIGURES/* .
  goto eof
)
if "%whichguides%" == "fdsvalg" (
  cd %svn_root%\Manuals\FDS_Validation_Guide\SCRIPT_FIGURES
  for /D %%d in (*) do (
      echo.
      echo copying files from %%d
      cd %%d

      Title Download FDS validation guide %%d images

      pscp %linux_logon%:%firebotrepo%/Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/%%d/* .
      cd ..
  )
  goto eof
)
if "%whichguides%" == "fdsverg" (
  Title Download FDS verification guide images

  cd %svn_root%\Manuals\FDS_Verification_Guide\SCRIPT_FIGURES
  pscp %linux_logon%:%firebotrepo%/Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/* .

  Title Download FDS verification guide scatterplot images

  cd %svn_root%\Manuals\FDS_Verification_Guide\SCRIPT_FIGURES\Scatterplots
  pscp %linux_logon%:%firebotrepo%/Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/Scatterplots/* .
  goto eof
)

:eof
echo.
echo copy complete
pause
