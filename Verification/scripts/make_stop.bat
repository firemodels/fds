@echo off

call %SVNROOT%\bot\Bundlebot\nightly\getopts.bat %*

echo 3 > %dir%\%infile%.stop
