@echo off

Rem Windows batch file to upload FDS release to the google download site.

Rem set version and revision to match those set in ZIP_FDS_Intel_win_64.bat

set version=5.3.0
set revision=3193

Rem uncomment 2 of the following 3 lines depending on what type of release is
Rem being uploaded

  Rem set level=Release-1_Major
  Rem set level=Release-2_Minor
  set level=Release-3_Maintenance

Rem directory containing googlecode password file (gc.passwd)
  set pwdir=c:\bin\

Rem -------- should not need to edit any lines below ---------


echo Uploading FDS %level% version=%version% revision=%revision%
pause

set upload_win64=1
Rem --------------- 64 bit Windows ----------------
if not %upload_win64% == 1 goto endif_win64
  set glabels=Type-Installer,Opsys-Windows,%level%
  set dplatform=64 bit Windows
  set platform=win64
  set summary=FDS Executable for %dplatform% (SVN r%revision%)
  set exe=fds_%version%_%platform%.exe
  echo -----
  echo Uploading %summary% - %exe%
  echo googlecode_upload.py --passwd-file-dir %pwdir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
       googlecode_upload.py --passwd-file-dir %pwdir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
:endif_win64

echo -----
echo Uploads complete
pause
