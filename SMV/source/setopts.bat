@echo off

echo arg1=%1
echo arg2=%2
set SIZE="-m32"
set COMPILER="gcc"
set COMPILER2="g++"

if "%1" NEQ "g" set COMPILER="icl"
if "%2" NEQ "3" set SIZE="-m64" 
