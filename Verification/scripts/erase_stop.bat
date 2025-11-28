@echo off

call %SVNROOT%\bot\Bundlebot\nightly\getopts.bat %*

if exist %dir%\%infile%.stop erase %dir%\%infile%.stop
