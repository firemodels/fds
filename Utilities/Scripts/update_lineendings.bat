@echo off
:: *** this script will not run if there are untracked or modified
::     files in your repository

:: This script uses commands found at:
:: https://help.github.com/articles/dealing-with-line-endings

set untracked_list=%temp%\untracked_list.txt
set untracked_count=%temp%\untracked_count.txt

set CURDIR=%CD%
cd ..\..
git clean -dxfn > %untracked_list% 
findstr /R /N "^.*" %untracked_list% | find /C ":" > %untracked_count%
set /p nuntracked=<%untracked_count%
erase %untracked_list% %untracked_count%
if %nuntracked% GTR 0 (
  echo.
  echo *** This repo has %nuntracked% untracked files.
  echo Clean the repo before proceeding.
  echo 1.  cd to the repo root 
  echo 2.  type: git clean -dxf
  goto eof
)

set dirtycount=%temp%\dirtycount.txt

git describe --long --dirty | find /C "dirty" > %dirtycount%
set /p ndirtycount=<%dirtycount%
erase %dirtycount%

if %ndirtycount% GTR 0 (
  echo.
  echo ***Warning: This repo has modified files.
  echo Commit or revert these changes before proceeding.
  echo Type: git status -uno
  echo to see which files have been changed.
  goto eof
fi

:: failsafe, should't get here if repo has untracked or modfied files
:: make sure repo is clean (otherwise untracked files will get committed)

echo Press enter to proceed with line ending update or CTRL c to abort
pause>Nul
git clean -dxf

git add . -u
git commit -m "saving files before refreshing line endings"
git rm --cached -r .
git reset --hard
git add .
git commit -m "normalize all the line endings"

:eof
cd %CURDIR%
