@echo off
set todir=%USERPROFILE%\bin
if NOT exist %todir% (
  echo creating %todir%
  echo add %todir% to your PATH variable
)
copy *.bat %todir%\.