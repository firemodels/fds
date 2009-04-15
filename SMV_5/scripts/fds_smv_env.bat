@echo off

Rem Configuration file for FDS-SMV scripts
Rem This file contains user customizable settings for the FDS and Smokeview build scripts


Rem VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
Rem VVVVVVVV SETTINGS BELOW CHANGED ONLY ONCE VVVVVVVVVVVVVVVVVVVVVVVVVVVV
Rem VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

Rem ---------- SVN repository locations ------------

set svn_root=d:\fds-smv
set svn_drive=d:

set linux_svn_root=FDS-SMV

Rem ---------- User/Host names -----------------

set linux_username=gforney

set google_username=gforney

set linux_hostname=acrux.cfr.nist.gov

Rem ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Rem ^^^^^^^^ SETTINGS ABOVE CHANGED ONLY ONCE ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Rem ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rem **********************************************************************

Rem ----------------------------------------------------------------------
Rem --- SETTINGS BELOW UPDATED MORE THAN ONCE - as repository changes ----
Rem ----------------------------------------------------------------------

Rem ----------- version and revision numbers -------------

set smv_version=test
set smv_revision=3758

set fds_version=5.3.1
set fds_revision=3729

Rem this parameter is not needed now but am keeping it around in case we change our minds
Rem about releasing docs on the download page

set docs_revision=3533

set verification_revision=3685

Rem ----------- for uploading to Google Code -----------------

Rem *** Smokeview - comment 2 of the following 3 lines
Rem    (the uncommented line is used for specifying type of Google upload)

Rem set smv_google_level=Release-1_Major
Rem set smv_google_level=Release-2_Minor
set     smv_google_level=Release-3_Maintenance

Rem *** Docs - comment 2 of the following 3 lines
Rem    (the uncommented line is used for specifying type of Google upload)

Rem set docs_google_level=Release-1_Major
Rem docs_google_level=Release-2_Minor
    set     docs_google_level=Release-3_Maintenance

Rem Example cases - comment 2 of the following 3 lines
Rem    (the uncommented line is used for specifying type of Google upload)

Rem set verification_google_level=Release-1_Major
Rem set verification_google_level=Release-2_Minor
set     verification_google_level=Release-3_Maintenance

Rem FDS - comment 2 of the following 3 lines
Rem    (the uncommented line is used for specifying type of Google upload)

Rem set fds_google_level=Release-1_Major
Rem set fds_google_level=Release-2_Minor
set     fds_google_level=Release-3_Maintenance


Rem --------------------------------------
Rem --- should not need to edit below ----
Rem --------------------------------------

set svn_logon=%linux_username%@%linux_hostname%

Rem put your google password into a file named gc.passwd located in
Rem   the directory below (your windows home directory found by simply
Rem   opening a winows command shell)

set google_password_dir=%homedrive%\%homepath%\

