@echo off
echo determine number of cmd.exe running in windows
echo after delaying 10 seconds
Timeout /t 10 >nul 
tasklist | find /i /c "cmd.exe" > temp.out
set /p numexe=<temp.out
echo numexe=%numexe%