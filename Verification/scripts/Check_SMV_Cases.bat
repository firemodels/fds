@echo off

set curdir=%CD%
set rungeomcases=1
set runwuicases=1
set runsmvcases=1

set SCRIPT_DIR=%CD%
cd ..
set BASEDIR="%CD%"
cd ..
set SVNROOT="%CD%"
cd %BASEDIR%

set stopscript=0
call :getopts %*
cd %curdir%
if %stopscript% == 1 (
  exit /b
)


set QFDS=call %SVNROOT%\Utilities\Scripts\checkfds.bat
set RUNCFAST=call %SVNROOT%\Utilities\Scripts\checkcfast.bat
set RUNTFDS=call %SVNROOT%\Utilities\Scripts\checkfds.bat

if "%runsmvcases%" == "1" (
  echo checking SMV cases
  call %SCRIPT_DIR%\SMV_Cases.bat
)
if "%rungeomcases%" == "1" (
  echo checking GEOM cases
  call %SCRIPT_DIR%\GEOM_Cases.bat
)
if "%runwuicases%" == "1" (
  echo checking WUI cases
  call %SCRIPT_DIR%\WUI_Cases.bat
)

goto eof

:getopts
 if (%1)==() exit /b
 set valid=0
 set arg=%1
 if /I "%1" EQU "-help" (
   call :usage
   set stopscript=1
   exit /b
 )
 if /I "%1" EQU "-geom" (
   set valid=1
   set runwuicases=0
   set runsmvcases=0
   set rungeomcases=1
 )
 if /I "%1" EQU "-smvwui" (
   set valid=1
   set runwuicases=1
   set runsmvcases=1
   set rungeomcases=0
 )
 if /I "%1" EQU "-wui" (
   set valid=1
   set runwuicases=1
   set runsmvcases=0
   set rungeomcases=0
 )
 shift
 if %valid% == 0 (
   echo.
   echo ***Error: the input argument %arg% is invalid
   echo.
   echo Usage:
   call :usage
   set stopscript=1
   exit /b
 )
if not (%1)==() goto getopts
exit /b

:usage  
echo Check_SMV_Cases [options]
echo. 
echo -help   - display this message
echo -geom   - run only geom cases
echo -smvwui - run only SMV and WUI cases
echo -wui    - run only WUI cases
exit /b


:eof
cd %curdir%

