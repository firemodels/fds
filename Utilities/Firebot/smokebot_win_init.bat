@echo   off

set CURDIR=%CD%

if NOT exist %userprofile%\smokebot (
  cd %userprofile%
  echo %userprofile%\smokebot does not exist - creating.
  svn co http://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Utilities/Firebot smokebot
  echo %userprofile%\smokebot created.
)

cd %CURDIR%
