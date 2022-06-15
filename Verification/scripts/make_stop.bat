@echo off

call %SVNROOT%\fds\Utilities\Scripts\getopts.bat %*

echo 3 > %dir%\%infile%.stop
