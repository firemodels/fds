@echo off

Rem Windows batch file to upload FDS release to the google download site.

Rem set version and revision to match those set in ZIP_FDS_Intel_win_32.bat

set version=upload_test
set revision=r3202

Rem uncomment 2 of the following 3 lines depending on what type of release is
Rem being uploaded

  Rem set level=Release-1_Major
  Rem set level=Release-2_Minor
  set level=Release-3_Maintenance

Rem directory containing googlecode password file (gc.passwd)
  set pwdir=d:\bin\

Rem -------- should not need to edit any lines below ---------


echo Uploading FDS %level% version=%version% revision=%revision%
pause

set upload_win32=1
Rem --------------- 32 bit Windows ----------------
if not %upload_win32% == 1 goto endif_win32
  set glabels=Type-Installer,Opsys-Windows,%level%
  set dplatform=32 bit Windows
  set platform=win32
  set summary=FDS %version% for %dplatform% (build %revision%)
  set exe=FDS_%version%_%revision%_%platform%.exe
  echo -----
  echo Uploading %summary% - %exe%
  echo googlecode_upload.py --passwd-file-dir %pwdir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
       googlecode_upload.py --passwd-file-dir %pwdir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
:endif_win32

echo -----
echo Uploads complete
pause
