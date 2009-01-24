@echo off
set version=test
set revision=3110

set release=1

if not %release% == 1 goto endif_release
  Rem set level=Release-1_Major
  Rem set level=Release-2_Minor
  set level=Release-3_Maintenance
  set upload_win32=1
  set upload_linux32=1
  set upload_osx32=1
:endif_release

if not %release% == 0 goto endif_test
  set level=Release-4_Test
  set upload_win32=1
  set upload_linux32=0
  set upload_osx32=0
:endif_test

echo Uploading Smokeview %level% version=%version% revision=%revision%
pause

cd ..\for_bundle\to_google
set pwdir=d:\bin\

Rem --------------- 32 bit Windows ----------------
if not %upload_win32% == 1 goto endif_win32
  set glabels=Type-Installer,Opsys-Windows,%level%
  set dplatform=32 bit Windows
  set platform=win32
  set summary=Smokeview %version% for %platform% build %revision%
  set exe=smv_%version%_%revision%_%platform%.exe
  echo -----
  echo Uploading %summary% - %exe%
Rem  echo googlecode_upload.py --passwd-file-dir %pwdir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
Rem       googlecode_upload.py --passwd-file-dir %pwdir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
:endif_win32

Rem --------------- 32 bit Linux ----------------

if not %upload_linux32% == 1 goto endif_linux32
  set glabels=Type-Archive,Opsys-Linux,%level%
  set dplatform=32 bit Linux
  set platform=linux
  set summary=Smokeview %version% for %dplatform% (build %revision%)
  set exe=smv_%version%_%revision%_%platform%.exe
  echo -----
  echo Uploading %summary% - %exe%
Rem  echo googlecode_upload.py --passwd-file-dir %pwdir% --config-dir none  -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
Rem       googlecode_upload.py --passwd-file-dir %pwdir% --config-dir none  -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
:endif_linux32

Rem --------------- 32 bit OSX ----------------

if not %upload_osx32% == 1 goto endif_osx32
  set glabels=Type-Archive,Opsys-OSX,%level%
  set dplatform=32 bit OSX
  set platform=osx
  set summary=Smokeview %version% for %dplatform% (build %revision%)
  set exe=smv_%version%_%revision%_%platform%.exe
  echo -----
  echo Uploading %summary% - %exe%
Rem  echo googlecode_upload.py --passwd-file-dir %pwdir% --config-dir none  -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
Rem       googlecode_upload.py --passwd-file-dir %pwdir% --config-dir none  -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
:endif_osx32

echo -----
echo Uploads complete
pause