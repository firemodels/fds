@echo   off

set CURDIR=%CD%

if NOT exist %userprofile%\smokebot (
  cd %userprofile%
  echo %userprofile%\smokebot does not exist - creating
  svn co http://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Utilities/Firebot smokebot
  echo %userprofile%\smokebot created.
)

:: create a clean cfast repository 

if NOT exist %userprofile%\cfastclean (
  cd %userprofile%
  echo %userprofile%\cfastclean does not exist - creating
  svn co http://cfast.googlecode.com/svn/trunk/cfast/trunk cfastclean
)

:: create a clean FDS repository

if NOT exist %userprofile%\FDS-SMVclean (
  cd %userprofile%
  echo %userprofile%\FDS-SMVclean does not exist - creating
  svn co http://fds-smv.googlecode.com/svn/trunk/FDS/trunk FDS-SMVclean
)

cd %CURDIR%
