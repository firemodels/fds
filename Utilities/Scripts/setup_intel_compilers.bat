@echo off

IF "%SETUP_IFORT_COMPILER_64%"=="1" GOTO envexist

  set SETUP_IFORT_COMPILER_64=1

  IF DEFINED IFORT_COMPILER14 set IFORT_COMPILER=%IFORT_COMPILER14%
  IF DEFINED IFORT_COMPILER15 set IFORT_COMPILER=%IFORT_COMPILER15%
  IF DEFINED IFORT_COMPILER16 set IFORT_COMPILER=%IFORT_COMPILER16%

  IF NOT DEFINED IFORT_COMPILER (
    echo "*** Error: Intel compiler environment variable not defined."
  )
  IF DEFINED IFORT_COMPILER (
    echo Setting up compiler environment
    call "%IFORT_COMPILER%\bin\compilervars" intel64
  )
:envexist
