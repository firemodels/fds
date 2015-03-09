@echo   off

set CURDIR=%CD%

if NOT exist %userprofile%\firebot (
  cd %userprofile%
  echo %userprofile%\firebot does not exist - creating.
  svn co http://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Utilities/Firebot firebot
  echo %userprofile%\firebot created.
)

cd %CURDIR%
