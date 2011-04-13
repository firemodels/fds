set intelbin="%IFORT_COMPILER11%\bin"

call %intelbin%\ifortvars ia32
if exist "%VS_COMPILER%\..\..\bin\vcvars32.bat" call "%VS_COMPILER%\..\..\bin\vcvars32.bat"

icl -o sh2bat sh2bat.c
pause
