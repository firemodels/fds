@echo   off

set CURDIR=%CD%

if NOT exist %userprofile%\firebot (
  cd %userprofile%
  echo %userprofile%\firebot does not exist - creating
  svn co http://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Utilities/Firebot firebot
  echo %userprofile%\firebot created.
)

:: create a clean FDS repository

if NOT exist %userprofile%\FDS-SMVclean (
  cd %userprofile%
  echo %userprofile%\FDS-SMVclean does not exist - creating
  svn co http://fds-smv.googlecode.com/svn/trunk/FDS/trunk FDS-SMVclean
)

cd %CURDIR%
