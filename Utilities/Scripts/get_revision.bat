@echo off
set svn_dir=%1

set svn_revision=Revision: unknown
set temp1=%temp%\temp.txt
set temp1c=%temp%\tempc.txt

set CURDIR=%CD%

:: looking for svn

svn 1> %temp1% 2>&1
type %temp1% | find /i /c "not recognized" > %temp1c%
set /p nothaveSVN=<%temp1c%
if %temp1c% == 1 (
  exit /b 1
)

:: looking for grep

grep 1> %temp1% 2>&1
type %temp1% | find /i /c "not recognized" > %temp1c%
set /p nothaveSVN=<%temp1c%
if %temp1c% == 1 (
  exit /b 1
)

cd %svn_dir%
svn info 2>&1 | grep Revision > %temp1%
set /p svn_revision=<%temp1%
cd %CURDIR%

