@echo off

:: ---- FDS and smokeview version strings for bundle
set fds_version=6.5.0
set smv_version=6.3.7

:: ---- FDS and smokeview revision version strings
::      (when creating test releases)
set smv_revision=3057a8c
set fds_revision=02a6097

:: ---- for obtaining log entries
set smvlogdate="01-Mar-2016"

:: PC repo

set svn_root=%userprofile%\FDS-SMV
set svn_drive=c:

:: Linux/OSX repo

set linux_svn_root=FDS-SMV
set firebotrepo=/home2/smokevis2/firebot/FDS-SMVgitclean
set smokebotrepo=/home2/smokevis2/smokebot/FDS-SMVgitclean
set fdswikirepo=%userprofile%\FDS-SMVwikis

:: FDS/Smokeview version

set fdssmv_major_version=6
set fds_edition=FDS6
set smv_edition=SMV6

:: Linux user and host name

set linux_hostname=blaze.nist.gov
set linux_username=%username%
set linux_logon=%linux_username%@%linux_hostname%

:: OSX user and host name

set osx_hostname=bluesky.el.nist.gov
set osx_username=%username%
set osx_logon=%osx_username%@%osx_hostname%
