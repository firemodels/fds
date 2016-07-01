@echo off

call %SVNROOT%\FDS\Utilities\Scripts\getopts.bat %*

if exist %dir%\%infile%.stop erase %dir%\%infile%.stop
