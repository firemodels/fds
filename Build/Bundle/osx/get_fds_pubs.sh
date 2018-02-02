#!/bin/bash

source ../scripts/GET_ENV.sh

if [ ! -e ~/FDS_Guides ]; then
  mkdir ~/FDS_Guides
fi

pubdir=$firebothome/.firebot/pubs
scp -q $linux_hostname\:$pubdir/FDS_Config_Management_Plan.pdf    ~/FDS_Guides/.
scp -q $linux_hostname\:$pubdir/FDS_Technical_Reference_Guide.pdf ~/FDS_Guides/.
scp -q $linux_hostname\:$pubdir/FDS_User_Guide.pdf                ~/FDS_Guides/.
scp -q $linux_hostname\:$pubdir/FDS_Validation_Guide.pdf          ~/FDS_Guides/.
scp -q $linux_hostname\:$pubdir/FDS_Verification_Guide.pdf        ~/FDS_Guides/.
