@echo off
set version=5.3.7
set revision=3110

Rem select level of release to be uploaded by adding/removing Rem's to/from 
Rem    the following lines

Rem set level=Release-1_Major
Rem set level=Release-2_Minor
set level=Release-3_Maintenance
Rem set level=Release-4_Test


Rem - do not edit below except to select which files to be uploaded 
Rem     (by adding/removing Rem's to/from googlecode_upload.py lines )

cd ..\for_bundle\to_google

set glabels=Type-Installer,Opsys-Windows,%level%
set platform=win32
set summary=smokeview %version% for %platform% (build %revision%)
set exe=smv_%version%_%revision%_%platform%.exe
echo googlecode_upload.py -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
rem remove rem on next line to upload to google
rem googlecode_upload.py -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%

set glabels=Type-Archive,Opsys-Linux,%level%
set platform=Linux
set summary=smokeview %version% for %platform% (build %revision%)
set exe=smv_%version%_%revision%_%platform%.exe
echo googlecode_upload.py -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
rem remove rem on next line to upload to google
rem googlecode_upload.py -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%


set glabels=Type-Archive,Opsys-OSX,%level%
set platform=OSX
set summary=smokeview %version% for %platform% (build %revision%)
set exe=smv_%version%_%revision%_%platform%.exe
echo googlecode_upload.py -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
rem remove rem on next line to upload to google
rem googlecode_upload.py -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
pause