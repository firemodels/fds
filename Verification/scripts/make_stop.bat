@echo off
call %SVNROOT%\Utilities\Scripts\getopts.bat %*

echo 5 > %dir%\%infile%.stop
