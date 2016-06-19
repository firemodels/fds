@echo off
:: *** warning: this script cleans your repo
::     DO NOT RUN if you have any uncommited work you wish to save

:: This script uses commands found at:
:: https://help.github.com/articles/dealing-with-line-endings

set untracked_list=untracked_list.txt
set untracked_count=untracked_count.txt

set CURDIR=%CD%
cd ..\..
git clean -dxfn > %untracked_list% 
findstr /R /N "^.*" %untracked_list% | find /C ":" > %untracked_count%
set /p nuntracked=<%untracked_count%
erase %untracked_list% %untracked_count%
if %nuntracked% GTR 0 (
  echo *** This repo has %nuntracked% untracked files.
  echo Before updating line endings, the repo must be cleaned.
  echo 1.  cd to the repo root 
  echo 2.  type: git clean -dxf
  goto eof
)

git clean -dxf
git add . -u
git commit -m "saving files before refreshing line endings"
git rm --cached -r .
git reset --hard
git add .
git commit -m "normalize all the line endings"

:eof
cd %CURDIR%
