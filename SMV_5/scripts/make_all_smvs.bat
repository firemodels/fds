@echo off
set version=5.3_2720
set logon=gforney@acrux.cfr.nist.gov

set scriptdir=FDS-SMV/SMV_5/scripts
set bundledir=FDS-SMV/SMV_5/for_bundle

plink %logon% %scriptdir%/svn_update.csh
plink %logon% %scriptdir%/make_smvs.csh
plink %logon% %scriptdir%/make_dists.csh %version%
pscp %logon%:%bundledir%/smv_%version%_osx.tar.gz ..\for_bundle\to_google\.
pscp %logon%:%bundledir%/smv_%version%_linux.tar.gz ..\for_bundle\to_google\.
pscp %logon%:%scriptdir%/make_linux.out ..\for_bundle\to_google\.
pscp %logon%:%scriptdir%/make_osx.out ..\for_bundle\to_google\.


make_smv_release_win32.bat %version%
pause
