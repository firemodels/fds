@echo off

Rem get FDS-SMV test cases for Windows

set REPOS=d:\fds-smv

Rem --------- should not need to edit below ----------

set revision=%1

cd %REPOS%\Utilities\Scripts\to_google\

set testdir=FDS_Test_cases_%revision%

svn export https://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Test_cases Test_cases
wzzip -a -r -P %testdir%.zip Test_cases
d:\bin\winzip\wzipse32 %testdir%.zip -d "c:\program files\nist\"
erase %testdir%.zip
pause
