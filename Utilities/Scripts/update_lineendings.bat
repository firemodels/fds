@echo off

set CURDIR=%CD%
cd ..\..
git clean -dxf
git add . -u
git commit -m "saving files before refreshing line endings"
git rm --cached -r .
git reset --hard
git add .
git commit -m "normalize all the line endings"
cd %CURDIR%
