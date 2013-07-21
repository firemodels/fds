@echo off
set SIZE=-m32
set SIZE2=ia32

set COMPILER=gcc
set COMPILER2=g++

set MSCOMPILER=X86

IF "%1" NEQ "g" SET COMPILER=icl
IF "%1" NEQ "g" SET COMPILER2=icl

IF "%1" EQU "m" SET COMPILER=cl
IF "%1" EQU "m" SET COMPILER2=cl

IF "%2" NEQ "3" SET SIZE=-m64
IF "%2" NEQ "3" SET SIZE2=intel64

IF "%2" NEQ "3" SET MSCOMPILER=x86_amd64

IF "%COMPILER%" == "icl" SET SIZE=
IF "%COMPILER%" == "cl" SET SIZE=

IF "%COMPILER%" NEQ "cl" GOTO MSenvexist

IF "%MSCOMPILERS_DEFINED%" EQU "1" GOTO MSenvexist
echo Setting up MS compiler environment
set MSCOMPILERS_DEFINED=1
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall" %MSCOMPILER%
GOTO Ienvexist
:MSenvexist

IF "%COMPILER%" NEQ "icl" GOTO Ienvexist
IF "%ICOMPILERS_DEFINED%" EQU "1" GOTO Ienvexist
echo Setting up Intel compiler environment
set ICOMPILERS_DEFINED=1
call "%IFORT_COMPILER13%\bin\compilervars" %SIZE2%
:Ienvexist
