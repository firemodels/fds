@echo off
Rem ------ FDS/Smokeview version and revision numbers ---------

set smv_version=v6test_8057
set smv_revision=8093
set fds_version=v6test_7978
set fds_revision=7978
set fds_edition=FDS5

Rem ---------- FDS-Smokeview repository settings ------------

set svn_root=%userprofile%\fds-smv
set svn_drive=c:
set linux_svn_root=FDS-SMV

Rem ---------- User/Host names -----------------

set linux_hostname=leo.cfr.nist.gov
set osx_hostname=bluesky.cfr.nist.gov
set LINUXCOMPILE=fire79

set linux_username=%username%
set svn_logon=%linux_username%@%linux_hostname%

Rem ----------- for uploading to Google Code -----------------

set google_username=%username%
set google_password_dir=%userprofile%\
set upload=%svn_root%\smv\scripts\googlecode_upload.py
set fds_google_level=Release-3_Maintenance
set smv_google_level=Release-3_Maintenance
