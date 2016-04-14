@echo off

:: Smokeview version and revision

set smv_version=6.3.4
set smv_revision=b22a887

:: FDS version and revision

set fds_version=test
set fds_revision=59edcfc

set smvlogdate="1-Nov-2015"

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

set osx_hostname=192.168.1.5
set osx_username=%username%
set osx_logon=%osx_username%@%osx_hostname%
