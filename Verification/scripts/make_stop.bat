@echo off
call %SVNROOT%\fds\Utilities\Scripts\getopts.bat %*

echo 2 > %dir%\%infile%.stop
