#!/bin/bash

source ../scripts/GET_ENV.sh

if [ ! -e ~/FDS_Guides ]; then
  mkdir ~/FDS_Guides
fi

pubdir=$smokebothome/.smokebot/pubs
scp $linux_hostname\:$pubdir/SMV_Technical_Reference_Guide.pdf ~/FDS_Guides/.
scp $linux_hostname\:$pubdir/SMV_User_Guide.pdf                ~/FDS_Guides/.
scp $linux_hostname\:$pubdir/SMV_Verification_Guide.pdf        ~/FDS_Guides/.
