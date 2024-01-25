@echo off

IF  X%SETVARS_COMPLETED% == X1 GOTO intel_envexist

  set "ONEAPIDIR=C:\Program Files (x86)\Intel\oneAPI"
  IF DEFINED ONEAPI_ROOT set "ONEAPIDIR=%ONEAPI_ROOT%"
  IF NOT EXIST "%ONEAPIDIR%\setvars.bat" goto intel_notexist

  echo Defining Intel compiler environment
  call "%ONEAPIDIR%\setvars" intel64
  set INTEL_IFORT=ifort

  IF  X%SETVARS_COMPLETED% == X1 GOTO intel_envexist

:intel_notexist
  echo ***error: Intel compiler environment is not setup
  goto :eof

:intel_envexist
:eof
