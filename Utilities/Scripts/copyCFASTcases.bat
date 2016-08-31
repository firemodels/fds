@echo off

call %SVNROOT%\fds\Utilities\Scripts\getopts.bat %*

:: set dir=%1
:: set infile=%2

set in=%infile%.in

if not exist %outdir%\%dir% mkdir %outdir%\%dir%
copy %dir%\%in% %outdir%\%dir%\.
if exist %dir%\%infile%.ini copy %dir%\%infile%.ini %outdir%\%dir%\.
if exist %dir%\%infile%.ssf copy %dir%\%infile%.ssf %outdir%\%dir%\.
