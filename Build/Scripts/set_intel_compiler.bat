@echo off
set arg1=%1
set arg2=%2

if NOT x%arg2% == x  goto endif1
if x%arg1%     == xbot   goto endif1
  set arg2=%arg1%
  set arg1=
:endif1

set INTEL_IFORT=%arg2%

if NOT x%arg2% == x goto endif2
   set INTEL_IFORT=ifort
:endif2

