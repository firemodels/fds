@echo off
if exist "%MKLROOT%\lib\intel64\mkl_blacs_intelmpi_lp64.lib" (
echo 1
exit /b
)
if exist "%MKLROOT%\lib\mkl_blacs_intelmpi_lp64.lib" (
echo 1
exit /b
)
echo 1