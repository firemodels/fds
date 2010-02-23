@echo off
Rem Configuration file for FDS-SMV scripts
Rem This file contains user customizable settings for the FDS and Smokeview build scripts

Rem VVVVVVVV SETTINGS BELOW CHANGED ONLY ONCE VVVVVVVVVVVVVVVVVVVVVVVVVVVV

Rem ---------- SVN repository locations ------------

set svn_root=d:\fds-smv
set svn_drive=d:
set linux_svn_root=FDS-SMV

Rem ---------- User/Host names -----------------

set linux_username=gforney
set google_username=gforney
set linux_hostname=acrux.cfr.nist.gov

Rem ^^^^^^^^ SETTINGS ABOVE CHANGED ONLY ONCE ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rem --- SETTINGS BELOW UPDATED MORE THAN ONCE - as repository changes ----

Rem ----------- version and revision numbers -------------

Rem ------ Smokeview version and revision numbers ---------

set smv_version=xyz
set smv_revision=5621

Rem ------ FDS version and revision numbers ---------

set fds_version=test
set fds_revision=5545

Rem ------ Verification case revision number ---------

set verification_revision=4584

Rem ------ including other documentation ---------

set docs_include_in_bundle=0

Rem ----------- for uploading to Google Code -----------------

Rem *** FDS - comment 2 of the following 3 lines
Rem    (the uncommented line is used for specifying type of Google upload)

Rem set fds_google_level=Release-1_Major
Rem set fds_google_level=Release-2_Minor
set     fds_google_level=Release-3_Maintenance

Rem *** Smokeview - comment 2 of the following 3 lines
Rem    (the uncommented line is used for specifying type of Google upload)

Rem set smv_google_level=Release-1_Major
set smv_google_level=Release-2_Minor
Rem set     smv_google_level=Release-3_Maintenance

Rem Example cases - comment 2 of the following 3 lines
Rem    (the uncommented line is used for specifying type of Google upload)

Rem set verification_google_level=Release-1_Major
Rem set verification_google_level=Release-2_Minor
set     verification_google_level=Release-3_Maintenance

Rem --------------------------------------
Rem --- should not need to edit below ----
Rem --------------------------------------

set svn_logon=%linux_username%@%linux_hostname%

set homedir=%homedrive%%homepath%\

Rem put your google password into a file named gc.passwd located in
Rem   the directory below (your windows home directory found by simply
Rem   opening a winows command shell)

set google_password_dir=d:\gpf_home\

