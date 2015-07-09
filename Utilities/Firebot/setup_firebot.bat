@echo   off

set CURDIR=%CD%
set gitrepo=FDS-SMVgitclean
set gitrepodir=%userprofile%\%gitrepo%
set botdir=%userprofile%\firebotgit

:: create a clean FDS repository

if NOT exist %gitrepodir% (
  cd %userprofile%
  echo %gitrepo% does not exist - creating
  git clone git@github.com:firemodels/fds-smv.git %gitrepo%
)

if NOT exist %botdir% (
  echo %botdir% does not exist - creating
  mkdir %botdir%
  cd %botdir%
  copy %gitrepodir%\Utilities\Firebot\*.bat
)

cd %CURDIR%
