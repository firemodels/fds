@echo off

Rem Configuration file for FDS-SMV scripts
Rem This file contains user customizable settings for the FDS and Smokeview build scripts


Rem ----------------------------------------------------------------------
Rem ------ settings below changed only once ------------------------------------
Rem ----------------------------------------------------------------------

Rem windows repository directory 
set svn_root=d:\fds-smv

Rem windows repository hard drive
set svn_drive=d:

Rem Linux repository directory 
set linux_svn_root=FDS-SMV

Rem Linux cluster user name
set linux_username=gforney

Rem google code user name
set google_username=gforney

Rem linux cluster hostname
set linux_hostname=acrux.cfr.nist.gov

Rem directory containing google password file
set google_password_dir=c:\bin\

Rem ----------------------------------------------------------------------
Rem ------ settings updated updated as repository changes ----------------
Rem ----------------------------------------------------------------------

Rem define Smokeview version and svn revision info

set smv_version=5.3.11
set smv_revision=3322

Rem define Google-code release level for Smokeview

Rem set smv_google_level=Release-1_Major
Rem set smv_google_level=Release-2_Minor
set     smv_google_level=Release-3_Maintenance

Rem define Google-code release level for docs

Rem set docs_google_level=Release-1_Major
Rem set docs_google_level=Release-2_Minor
set     docs_google_level=Release-3_Maintenance

Rem define Google-code release level for Test_cases

Rem set test_cases_google_level=Release-1_Major
Rem set test_cases_google_level=Release-2_Minor
set     test_cases_google_level=Release-3_Maintenance

Rem define FDS version and svn revision info

set fds_version=uploadtest
set fds_revision=donotuse

Rem define Google-code release level for FDS

Rem set fds_google_level=Release-1_Major
Rem set fds_google_level=Release-2_Minor
set     fds_google_level=Release-3_Maintenance


Rem revision number for Documentation directory

set docs_revision=3235

Rem define revision number for test cases

set test_cases_revision=3222a

Rem --------------------------------------
Rem --- do not edit below ----------------
Rem --------------------------------------

set svn_logon=%linux_username%@%linux_hostname%
