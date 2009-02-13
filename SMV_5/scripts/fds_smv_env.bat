@echo off

Rem Configuration file for FDS-SMV scripts
Rem This file contains user customizable settings for the FDS and Smokeview build scripts


Rem ----------------------------------------------------------------------
Rem ------ settings below changed only once ------------------------------------
Rem ----------------------------------------------------------------------

Rem repository directory 
set svn_root=d:\fds-smv

Rem repository hard drive
set svn_drive=d:

Rem Linux cluster user name
set linux_username=gforney

Rem google code user name
set google_username=gforney

Rem linux cluster hostname
set linux_hostname=acrux.cfr.nist.gov

Rem ----------------------------------------------------------------------
Rem ------ settings updated updated as repository changes ----------------
Rem ----------------------------------------------------------------------

Rem define Smokeview version and svn revision info

set smv_version=5.3.11
set smv_revision=3279

Rem define FDS version and svn revision info

set fds_version=xxxx
set fds_revision=zzzz

Rem revision number for Documentation directory

set docs_revision=3235

Rem define revision number for test cases

set test_cases_revision=3222a

Rem --------------------------------------
Rem --- do not edit below ----------------
Rem --------------------------------------

set svn_logon=%linux_username%@%linux_hostname%
