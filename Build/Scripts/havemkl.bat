@echo off
if not exist "%MKLROOT%\lib\intel64\mkl_blacs_intelmpi_lp64.lib" (
echo 0
exit /b
)
echo 1