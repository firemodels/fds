@echo off

set curdir=%CD%
set running=bot.running
if not exist %running% (
  svn update
  echo 1 > %running%
  call smokebot_win_git.bat 1
  cd %curdir%
  erase %running%
) else (
  echo A bot is already running.
  echo Erase the file %running% if this is not the case
)
