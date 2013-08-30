@echo off
Rem ------ FDS/Smokeview version and revision numbers ---------

set smv_version=6.0.6
set smv_revision=13408

set fds_version=bundletest9
set fds_revision=13139

set fdssmv_major_version=6

Rem ---------- App names ------------

set wzzip=wzzip

Rem ---------- FDS-Smokeview repository settings ------------

set fds_edition=FDS6
set svn_root=%userprofile%\fds-smv
set svn_drive=c:
set linux_svn_root=FDS-SMV

Rem ---------- User/Host names -----------------

set linux_hostname=blaze.nist.gov
set osx_hostname=bluesky

set linux_username=%username%
set svn_logon=%linux_username%@%linux_hostname%

Rem ----------- for uploading to Google Code -----------------

set google_username=%username%
set google_password_dir=%userprofile%\
set upload=%svn_root%\smv\scripts\googlecode_upload.py
set fds_google_level=Release-3_Maintenance
set smv_google_level=Release-3_Maintenance
