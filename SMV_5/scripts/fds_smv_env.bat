@echo off

Rem Configuration file for FDS-SMV scripts
Rem This file contains user customizable settings for the FDS and Smokeview build scripts


Rem ----------------------------------------------------------------------
Rem ------ settings below changed only once ------------------------------------
Rem ----------------------------------------------------------------------

Rem ---------- SVN repository locations ------------

set svn_root=d:\fds-smv
set svn_drive=d:

set linux_svn_root=FDS-SMV

Rem ---------- User/Host names -----------------

set linux_username=gforney

set google_username=gforney

set linux_hostname=acrux.cfr.nist.gov

Rem ----------------------------------------------------------------------
Rem ------ settings updated updated as repository changes ----------------
Rem ----------------------------------------------------------------------

Rem ----------- Software version and revision numbers -------------

set smv_version=5.3.11
set smv_revision=3322

set fds_version=xxxx
set fds_revision=yyyy

set docs_revision=3371

set verification_revision=3222a


Rem ----------- For Google Code Download site -----------------

Rem *** Smokeview - comment 2 of the following 3 lines

Rem set smv_google_level=Release-1_Major
Rem set smv_google_level=Release-2_Minor
set     smv_google_level=Release-3_Maintenance

Rem *** Docs - comment 2 of the following 3 lines

Rem set docs_google_level=Release-1_Major
Rem docs_google_level=Release-2_Minor
    set     docs_google_level=Release-3_Maintenance

Rem Example cases - comment 2 of the following 3 lines

Rem set verification_google_level=Release-1_Major
Rem set verification_google_level=Release-2_Minor
set     verification_google_level=Release-3_Maintenance

Rem define Google-code release level for FDS
Rem  *** comment 2 of the following 3 lines

Rem set fds_google_level=Release-1_Major
Rem set fds_google_level=Release-2_Minor
set     fds_google_level=Release-3_Maintenance


Rem --------------------------------------
Rem --- should not need to edit below ----------------
Rem --------------------------------------

set svn_logon=%linux_username%@%linux_hostname%

Rem put your google password into a file named gc.passwd located in
Rem   the directory below (your windows home directory found by simply
Rem   opening a winows command shell)

set google_password_dir=%homedrive%\%homepath%\

