@echo off
set ARG=%1

set INTEL_IFORT=%ARG%
if NOT x%ARG% == x goto endif1
   set INTEL_IFORT=ifort
:endif1

