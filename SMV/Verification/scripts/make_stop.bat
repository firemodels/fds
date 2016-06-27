@echo off
call %SVNROOT%\Utilities\Scripts\getopts.bat %*

echo 2 > %dir%\%infile%.stop
