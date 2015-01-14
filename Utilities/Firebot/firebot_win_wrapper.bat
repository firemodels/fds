@echo off
if not exist firebot_win_running.txt (
  echo 1 > firebot_win_running.txt
  svn update
  call firebot_win.bat
  erase firebot_win_running.txt
)
