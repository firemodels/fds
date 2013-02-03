@echo off

echo arg1=%1
echo arg2=%2
set SIZE="-m32"
set COMPILER="gcc"
set COMPILER2="g++"

if "%1" NEQ "g" set COMPILER="icc"
if "%1" NEQ "g" set COMPILER2="icc"
if "%2" NEQ "3" set SIZE="-m64" 
