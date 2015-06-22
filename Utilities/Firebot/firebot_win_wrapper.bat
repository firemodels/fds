@echo off

set emailto=%1

set curdir=%CD%
set running=bot.running
if not exist %running% (
  svn update
  echo 1 > %running%
  call firebot_win.bat %emailto%
  cd %curdir%
  erase %running%
) else (
  echo A bot is already running.
  echo Erase the file %running% if this is not the case
)
