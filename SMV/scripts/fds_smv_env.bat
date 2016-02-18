@echo off
:: ------ FDS/Smokeview version and revision numbers ---------

set smv_version=6.3.4
set smv_revision=5c664c0

set fds_version=test
set fds_revision=59edcfc

set fdssmv_major_version=6
set smvdate="1-Nov-2015"
set githash=5c2cfb1

:: ---------- FDS-Smokeview repository settings ------------

set fds_edition=FDS6
set smv_edition=SMV6
set svn_root=%userprofile%\FDS-SMV
set svn_drive=c:
set linux_svn_root=FDS-SMV
set firebotrepo=/home2/smokevis2/firebot/FDS-SMVgitclean
set smokebotrepo=/home2/smokevis2/smokebot/FDS-SMVgitclean
set fdswikirepo=%userprofile%\FDS-SMVwikis

:: ---------- User/Host names -----------------

set linux_hostname=blaze.nist.gov
set linux_username=%username%
set linux_logon=%linux_username%@%linux_hostname%

set osx_hostname=bluesky.el.nist.gov
set osx_username=%username%
set osx_logon=%osx_username%@%osx_hostname%

set svn_logon=%linux_username%@%linux_hostname%
