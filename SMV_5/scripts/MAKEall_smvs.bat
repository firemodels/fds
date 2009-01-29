@echo off

Rem  Windows batch file to build Smokeview for all platforms.
Rem  This script builds LInux and OSX Smokeview's by doing a
Rem  remote shell (plink) to the NIST Linux cluster. 

set version=5.3.9_3167
Rem set version=test

Rem -----------------------------------------------------------
Rem shouldn't need to change any lines below

set logon=gforney@acrux.cfr.nist.gov

set scriptdir=FDS-SMV/SMV_5/scripts
set bundledir=FDS-SMV/SMV_5/for_bundle

plink %logon% %scriptdir%/svn_update.csh
plink %logon% %scriptdir%/make_smvs.csh
plink %logon% %scriptdir%/make_dists.csh %version%

echo downloading Linux Smokeview files
pscp %logon%:%bundledir%/smv_%version%_linux.tar.gz ..\for_bundle\to_google\.
pscp %logon%:%bundledir%/smv_%version%_linux_64.tar.gz ..\for_bundle\to_google\.
pscp %logon%:%scriptdir%/make_intel_linux_32.out ..\for_bundle\to_google\.
pscp %logon%:%scriptdir%/make_intel_linux_64.out ..\for_bundle\to_google\.

echo downloading MAC OSX Smokeview files
pscp %logon%:%bundledir%/smv_%version%_osx.tar.gz ..\for_bundle\to_google\.
pscp %logon%:%scriptdir%/make_intel_osx_32.out ..\for_bundle\to_google\.

call make_smv_release_win32.bat %version%

echo uploading Windows Smokeview files
pscp  ..\for_bundle\smokeview_release.exe %logon%:%bundledir%/.
pscp  ..\for_bundle\smokezip_release.exe %logon%:%bundledir%/.
pause
