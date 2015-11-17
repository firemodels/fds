@echo off
:: for now only support -d and one non-option parameter (filename)

:GETOPTS
 if /I "%1" EQU "-d" (
   set dir=%2
   shift
 )
 if /I "%1" EQU "-p" (
   set nprocs=%2
   shift
 )
 set infile=%1
 shift
if not (%1)==() goto GETOPTS

:: strip extension
for %%f in (%infile%) do (
    set infile=%%~nf
)
