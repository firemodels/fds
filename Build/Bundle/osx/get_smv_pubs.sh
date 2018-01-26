#!/bin/bash

if [ -e ~/FDS_SMV_ENVpc.sh ]; then
  source ~/FDS_SMV_ENVpc.sh
else
  source ~/FDS_SMV_ENV.sh
fi

if [ ! -e ~/FDS_Guides ]; then
  mkdir ~/FDS_Guides
fi

pubdir=$smokebothome/.smokebot/pubs
scp -q $linux_hostname\:$pubdir/SMV_Technical_Reference_Guide.pdf ~/FDS_Guides/.
scp -q $linux_hostname\:$pubdir/SMV_User_Guide.pdf                ~/FDS_Guides/.
scp -q $linux_hostname\:$pubdir/SMV_Verification_Guide.pdf        ~/FDS_Guides/.
