@echo off

call %SVNROOT%\Utilities\Scripts\getopts.bat %*

if exist %dir%\%infile%.stop erase %dir%\%infile%.stop