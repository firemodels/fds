@echo off
set curdir=%CD%
set running=bot.running
if not exist %running% (
  echo 1 > %running%
  call smokebot_win.bat 0
  cd %curdir%
  erase %running%
) else (
  echo A bot is already running.
  echo Erase the file %running% if this is not the case
)
