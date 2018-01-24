#!/bin/bash

if [ -e ~/FDS_SMV_ENVpc.sh ]; then
  source ~/FDS_SMV_ENVpc.sh
else
  source ~/FDS_SMV_ENV.sh
fi

if [ ! -e ~/FDS_Guides ]; then
  mkdir ~/FDS_Guides
fi

pubdir=$firebotrepo/fds/Manuals
scp $linux_hostname\:$pubdir/FDS_Config_Management_Plan/FDS_Config_Management_Plan.pdf       ~/FDS_Guides/.
scp $linux_hostname\:$pubdir/FDS_Technical_Reference_Guide/FDS_Technical_Reference_Guide.pdf ~/FDS_Guides/.
scp $linux_hostname\:$pubdir/FDS_User_Reference_Guide/FDS_User_Reference_Guide.pdf           ~/FDS_Guides/.
scp $linux_hostname\:$pubdir/FDS_Validation_Guide/FDS_Validation_Guide.pdf                   ~/FDS_Guides/.
scp $linux_hostname\:$pubdir/FDS_Verification_Guide/FDS_Verification_Guide.pdf               ~/FDS_Guides/.
