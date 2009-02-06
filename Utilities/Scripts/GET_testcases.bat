@echo off

Rem get FDS-SMV test cases for current revision

set revision=3222
set REPOS=d:\fds-smv

cd %REPOS%\Utilities\Scripts\to_google\

set testdir=FDS_Test_cases_%revision%

svn export https://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Test_cases %testdir%
cd %testdir%
wzzip -a -r -P %testdir%.zip *
d:\bin\winzip\wzipse32 %testdir%.zip -d "c:\program files\nist\"
copy %testdir%.exe ..\.
erase %testdir%.zip
erase %testdir%.exe
pause
