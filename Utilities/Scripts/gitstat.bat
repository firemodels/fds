@echo off

set gitout=gitout.txt
set gitcount=gitcount.txt

:: make sure git is installed

git 1> %gitout% 2>&1
type %gitout% | find /i /c "not recognized" > %gitcount%
set /p nothaveGIT=<%gitcount%
if %nothaveGIT% == 1 (
  echo 
  erase %gitout% %gitcount%
  exit /b 1
)
erase %gitout% %gitcount%

git diff --shortstat ..\..\Source 1> %gitout% 2>&1
type %gitout% | find /i /c "changed" > %gitcount%
set /p changes=<%gitcount%
if %changes% == 1 (
  echo -dirty
  erase %gitout% %gitcount%
  exit /b 1
)
erase %gitout% %gitcount%
