@echo off

:: $Date$ 
:: $Revision$
:: $Author$

call %SVNROOT%\fds\Utilities\Scripts\getopts.bat %*

set fulldir=%BASEDIR%/%dir%

set errfile=%infile%.err
set err2file=%infile%.err2
set err3file=%infile%.err3
set nerrorfile=%infile%.nerr

cd %fulldir%

grep -v "commands for target" %errfile% > %OUTDIR%\stage_error0.txt
grep -v "mpif.h" %OUTDIR%\stage_error0.txt > %err2file%
grep -H -i -A 5 -B 5 "error"           %err2file% >  %err3file%
grep -H -i -A 5 -B 5 "forrtl"          %err2file% >> %err3file%
grep -H -i -A 5 -B 5 "Run aborted"     %err2file% >> %err3file%
grep -H -i -A 5 -B 5 "STOP: Numerical" %err2file% >> %err3file%
grep -H -i -A 5 -B 5 "Segmentation"    %err2file% >> %err3file%
type %err3file% | find /v /c "  "> %nerrorfile%
set /p nerrors=<%nerrorfile%
if %nerrors% GTR 0 (
   type %err3file% >> %OUTDIR%\stage_error.txt
   echo. >> %OUTDIR%\stage_error.txt
)

cd %BASEDIR%
