@echo off
Rem ------ FDS/Smokeview version and revision numbers ---------

set smv_version=test
set smv_revision=59edcfc

set fds_version=test
set fds_revision=59edcfc

set fdssmv_major_version=6

Rem ---------- FDS-Smokeview repository settings ------------

set fds_edition=FDS6
set smv_edtion=SMV6
set svn_root=%userprofile%\fds-smv
set svn_drive=c:
set linux_svn_root=FDS-SMV

Rem ---------- User/Host names -----------------

set linux_hostname=blaze.nist.gov
set osx_hostname=bluesky.el.nist.gov

set linux_username=%username%
set svn_logon=%linux_username%@%linux_hostname%

Rem ----------- for uploading to Bintray -----------------

set bintray_api_key=%userprofile%\keys\bintray_api_key.txt
set upload=%svn_root%\Utilities\Scripts\curl
