@echo off

:: $Date$ 
:: $Revision$
:: $Author$

set curdir=%CD%
set running=firebot_win_running.status
if not exist %running% (
  svn update
  echo 1 > %running%
  call firebot_win.bat debug
  cd %curdir%
  erase %running%
) else (
  echo firebot_win is already running
  echo erase the file %running% if this is not the case
)
