@echo off

:: $Date$ 
:: $Revision$
:: $Author$

set curdir=%CD%
set running=smokebot_win_running.status
if not exist %running% (
  svn update
  echo 1 > %running%
  call smokebot_win.bat debug
  cd %curdir%
  erase %running%
) else (
  echo smokebot_win is already running
  echo erase the file %running% if this is not the case
)
