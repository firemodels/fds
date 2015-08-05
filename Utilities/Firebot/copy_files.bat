@echo off
set CURDIR=%CD%
set repo=FDS-SMVgitclean

cd %userprofile%\%repo%
git remote update
git pull
cd %CURDIR%
copy %userprofile%\%repo%\Utilities\Firebot\*.bat %CURDIR%
