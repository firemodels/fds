@echo off

IF "%SETUP_IFORT_COMPILER_64%"=="1" GOTO envexist

  set SETUP_IFORT_COMPILER_64=1

  IF DEFINED IFORT_COMPILER14 set IFORT_COMPILER=%IFORT_COMPILER14%
  IF DEFINED IFORT_COMPILER15 set IFORT_COMPILER=%IFORT_COMPILER15%
  IF DEFINED IFORT_COMPILER16 set IFORT_COMPILER=%IFORT_COMPILER16%
  IF DEFINED IFORT_COMPILER17 set IFORT_COMPILER=%IFORT_COMPILER17%
  IF DEFINED IFORT_COMPILER18 set IFORT_COMPILER=%IFORT_COMPILER18%
  IF DEFINED IFORT_COMPILER19 set IFORT_COMPILER=%IFORT_COMPILER19%
  IF DEFINED IFORT_COMPILER20 set IFORT_COMPILER=%IFORT_COMPILER20%

  IF NOT DEFINED IFORT_COMPILER (
    echo "*** Error: Intel compiler environment variable, IFORT_COMPILER, not defined."
    echo "    Intel compilers probably not installed."
    exit /b
  )
  IF NOT DEFINED I_MPI_ROOT (
    echo "*** Error: Intel MPI environment variable, I_MPI_ROOT, not defined."
    echo "    Intel MPI development environment probably not installed."
    exit /b
  )

  echo Setting up compiler environment
  set STARTUP="%IFORT_COMPILER%\bin\compilervars"
  set "I_MPI_ROOT_SAVE=%I_MPI_ROOT%"
  call %STARTUP% intel64

  echo Setting up MPI environment
  set "I_MPI_ROOT=%I_MPI_ROOT_SAVE%"
  set I_MPI_RELEASE_ROOT=%I_MPI_ROOT%\intel64\lib
  set   I_MPI_DEBUG_ROOT=%I_MPI_ROOT%\intel64\lib
  IF DEFINED IFORT_COMPILER19 set I_MPI_RELEASE_ROOT=%I_MPI_ROOT%\intel64\lib\release
  IF DEFINED IFORT_COMPILER19 set I_MPI_DEBUG_ROOT=%I_MPI_ROOT%\intel64\lib\debug
  call "%I_MPI_ROOT%\intel64\bin\mpivars" release

:envexist
