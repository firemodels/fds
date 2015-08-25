@echo off
git 1> git.out 2>&1
type git.out | find /i /c "not recognized" > git_count.txt
set /p nothaveGIT=<git_count.txt
if %nothaveGIT% == 1 (
  echo unknown
  erase git.out git_count.txt
  exit /b 1
)
erase git.out git_count.txt
git log -1 --format=%%cd