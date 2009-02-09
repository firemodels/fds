@echo off

Rem Configuration file for FDS-SMV scripts
Rem This file contains user customizable settings for the FDS and Smokeview build scripts


Rem ----------------------------------------------------------------------
Rem ------ settings updated updated as repository changes ----------------
Rem ----------------------------------------------------------------------

Rem define Smokeview version and svn revision info

set smv_version=zzzz
set smv_revision=qqqq

Rem define FDS version and svn revision info

set fds_version=xxxx
set fds_revision=zzzz

Rem revision number for Documentation directory

set docs_revision=wwww

Rem define revision number for test cases

set test_cases_revision=tttt

Rem ----------------------------------------------------------------------
Rem ------ settings changed only once ------------------------------------
Rem ----------------------------------------------------------------------

Rem define svn_root and svn_drive to point to FDS-SMV svn repository

set svn_root=d:\fds-smv
set svn_drive=d:

Rem set username for Linux cluster and google code

set linux_username=username
set google_username=username

Rem set linux cluster hostname

set linux_hostname=hostname

Rem --------------------------------------
Rem --- do not edit below ----------------
Rem --------------------------------------

set svn_logon=%linux_username%@%linux_hostname%
