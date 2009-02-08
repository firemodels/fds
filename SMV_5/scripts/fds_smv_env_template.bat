@echo off

Rem Define the environment variables below and copy this file to c:\bin renaming it fds_smv_env.bat


Rem define svn_root and svn_drive to point to FDS-SMV svn repository

set svn_root=d:\fds-smv
set svn_drive=d:

Rem set username for Linux cluster and google code

set linux_username=gforney
set google_username=gforney

Rem set linux cluster hostname

set linux_hostname=acrux.cfr.nist.gov

Rem define Smokeview version and svn revision info

set smv_version=ssss
set smv_revision=tttt

Rem define FDS version and svn revision info

set fds_version=xxxx
set fds_revision=zzzz


set svn_logon=%linux_username%@%linux_hostname%
