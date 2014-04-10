@echo off

:: example of using get_time and subtracting two dos variables

get_time_64 > time.txt
set /p time1=<time.txt
set base=100000000
set /a diff=%time1%-%base%
echo time1=%time1%
echo base=%base%
echo diff=%diff%