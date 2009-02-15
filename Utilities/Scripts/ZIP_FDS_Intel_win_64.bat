@echo off

Rem Batch file used to create a self-extracting archive containing FDS

set version=5.3.0
set revision=3193

set bin=c:\bin

Rem --------- should not need to edit lines below ------

set togoogle=to_google
set fdsdir=fds_%version%_win64

copy ..\Makefile\Intel_Win_64\fds5_win_64.exe %togoogle%\fds5_win64.exe

echo -
echo winzipping distribution directory
cd %togoogle%
if EXIST %fdsdir%.zip erase %fdsdir%.zip
%bin%\wzzip -a -r -P %fdsdir%.zip fds5_win64.exe

echo -
echo creating self-extracting archive
if EXIST %fdsdir%.exe erase %fdsdir%.exe
%bin%\winzip\wzipse32 %fdsdir%.zip -d "C:\Program Files (x86)\nist\FDS"
copy %fdsdir%.exe ..\.

pause