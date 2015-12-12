@echo off
Rem ------ FDS/Smokeview version and revision numbers ---------

set smv_version=6.3.2
set smv_revision=3988zc8

set fds_version=6.3.1
set fds_revision=59edcfc

set fdssmv_major_version=6
set smvdate="1-Nov-2015"
set githash=5c2cfb1

Rem ---------- FDS-Smokeview repository settings ------------

set fds_edition=FDS6
set smv_edition=SMV6
set svn_root=%userprofile%\FDS-SMV
set svn_drive=c:
set linux_svn_root=FDS-SMV
set firebotrepo=/home2/smokevis2/firebot/FDS-SMVgitclean
set smokebotrepo=/home2/smokevis2/smokebot/FDS-SMVgitclean

Rem ---------- User/Host names -----------------

set linux_hostname=blaze.nist.gov
set osx_hostname=floga.el.nist.gov

set linux_username=%username%
set svn_logon=%linux_username%@%linux_hostname%

Rem ----------- for uploading to Bintray -----------------

set bintray_api_key=%userprofile%\keys\bintray_api_key.txt
set upload=%svn_root%\Utilities\Scripts\curl
