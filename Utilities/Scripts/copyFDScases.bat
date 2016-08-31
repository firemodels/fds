@echo off

call %SVNROOT%\fds\Utilities\Scripts\getopts.bat %*

set in=%infile%.fds

if not exist %outdir%\%dir% mkdir %outdir%\%dir%
copy %dir%\%in% %outdir%\%dir%\.
if not exist %dir%\%in% echo %in%
if exist %dir%\%infile%.ini copy %dir%\%infile%.ini %outdir%\%dir%\.
if exist %dir%\%infile%.ssf copy %dir%\%infile%.ssf %outdir%\%dir%\.
