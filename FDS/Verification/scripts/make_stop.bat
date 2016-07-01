@echo off
call %SVNROOT%\FDS\Utilities\Scripts\getopts.bat %*

echo 2 > %dir%\%infile%.stop
