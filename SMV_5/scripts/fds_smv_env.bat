@echo off
Rem ---------- FDS-Smokeview repository settings ------------

set svn_root=d:\fds-smv
set svn_drive=d:
set linux_svn_root=FDS-SMV


Rem ------ FDS/Smokeview version and revision numbers ---------

set smv_version=5.5
set smv_revision=5868
set fds_version=5.5
set fds_revision=5545

Rem ---------- User/Host names -----------------

set linux_username=%username%
set google_username=%username%
set linux_hostname=acrux.cfr.nist.gov
set OSXHOST=tiger.cfr.nist.gov

Rem ----------- for uploading to Google Code -----------------

Rem google password file
set google_password_dir=d:\gpf_home\

set fds_google_level=Release-3_Maintenance
set smv_google_level=Release-3_Maintenance

Rem linux logon

set svn_logon=%linux_username%@%linux_hostname%




