set intelbin="%IFORT_COMPILER11%\bin"

call %intelbin%\ifortvars ia32

icl -o sh2bat sh2bat.c
pause
