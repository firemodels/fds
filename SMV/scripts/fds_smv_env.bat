@echo off
Rem ------ FDS/Smokeview version and revision numbers ---------

set smv_version=6.0.1
set smv_revision=11112

set fds_version=6.0bundletest
set fds_revision=11109

set fdssmv_major_version=6

Rem ---------- FDS-Smokeview repository settings ------------

set fds_edition=FDS%fdssmv_major_version%
set svn_root=%userprofile%\fds-smv
set svn_drive=c:
set linux_svn_root=FDS-SMV

Rem ---------- User/Host names -----------------

set linux_hostname=blaze.nist.gov
set osx_hostname=bluesky.cfr.nist.gov

set linux_username=%username%
set svn_logon=%linux_username%@%linux_hostname%

Rem ----------- for uploading to Google Code -----------------

set google_username=%username%
set google_password_dir=%userprofile%\
set upload=%svn_root%\smv\scripts\googlecode_upload.py
set fds_google_level=Release-3_Maintenance
set smv_google_level=Release-3_Maintenance
