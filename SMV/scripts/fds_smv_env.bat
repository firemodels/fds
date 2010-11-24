c@echo off
Rem ------ FDS/Smokeview version and revision numbers ---------

set smv_version=5.6
set smv_revision=7032
set fds_version=5.5.3
set fds_revision=7031

Rem ---------- FDS-Smokeview repository settings ------------

set svn_root=d:\fds-smv
set svn_drive=d:
set linux_svn_root=FDS-SMV

Rem ---------- User/Host names -----------------

set linux_hostname=acrux.cfr.nist.gov
set OSXHOST=bluesky.cfr.nist.gov
set LINUXCOMPILE=fire79

set linux_username=%username%
set google_username=%username%
set svn_logon=%linux_username%@%linux_hostname%

Rem ----------- for uploading to Google Code -----------------

set google_password_dir=%HOMEDRIVE%%HOMEPATH%\
set upload==%svn_root%\smv\scripts\googlecode_upload.py
set fds_google_level=Release-3_Maintenance
set smv_google_level=Release-3_Maintenance
