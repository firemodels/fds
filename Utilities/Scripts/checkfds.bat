@echo off

:: $Date$ 
:: $Revision$
:: $Author$

call %SVNROOT%\Utilities\Scripts\getopts.bat %*

set fulldir=%BASEDIR%/%dir%

set errfile=%infile%.err
set err2file=%infile%.err2

cd %fulldir%

grep -v "commands for target" %errfile% > %OUTDIR%\stage_error0.txt
grep -v "mpif.h" %OUTDIR%\stage_error0.txt > %err2file%
grep -H -i -A 5 -B 5 "error"           %err2file% >> %OUTDIR%\stage_error.txt
grep -H -i -A 5 -B 5 "forrtl"          %err2file% >> %OUTDIR%\stage_error.txt
grep -H -i -A 5 -B 5 "Run aborted"     %err2file% >> %OUTDIR%\stage_error.txt
grep -H -i -A 5 -B 5 "STOP: Numerical" %err2file% >> %OUTDIR%\stage_error.txt
grep -H -i -A 5 -B 5 "Segmentation"    %err2file% >> %OUTDIR%\stage_error.txt
echo. >> %OUTDIR%\stage_error.txt

cd %BASEDIR%
