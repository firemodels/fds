@echo off
set version=5.3.8_3138
set logon=gforney@acrux.cfr.nist.gov

set scriptdir=FDS-SMV/SMV_5/scripts
set bundledir=FDS-SMV/SMV_5/for_bundle

plink %logon% %scriptdir%/svn_update.csh
plink %logon% %scriptdir%/make_smvs.csh
plink %logon% %scriptdir%/make_dists.csh %version%

echo downloading Linux and OSX Smokeview files

pscp %logon%:%bundledir%/smv_%version%_osx.tar.gz ..\for_bundle\to_google\.
pscp %logon%:%bundledir%/smv_%version%_linux.tar.gz ..\for_bundle\to_google\.
pscp %logon%:%bundledir%/smv_%version%_linux_64.tar.gz ..\for_bundle\to_google\.
pscp %logon%:%scriptdir%/make_intel_linux_32.out ..\for_bundle\to_google\.
pscp %logon%:%scriptdir%/make_intel_linux_64.out ..\for_bundle\to_google\.
pscp %logon%:%scriptdir%/make_intel_osx_32.out ..\for_bundle\to_google\.

call make_smv_release_win32.bat %version%

echo uploading Windows Smokeview files
pscp  ..\for_bundle\smokeview_release.exe %logon%:%bundledir%/.
pscp  ..\for_bundle\smokezip_release.exe %logon%:%bundledir%/.
pause
