#!/bin/bash

source ../scripts/GET_ENV.sh

mkdir -p $GUIDE_DIR
pubdir=$smokebothome/.smokebot/pubs
scp -q $linux_hostname\:$pubdir/SMV_Technical_Reference_Guide.pdf $GUIDE_DIR/.
scp -q $linux_hostname\:$pubdir/SMV_User_Guide.pdf                $GUIDE_DIR/.
scp -q $linux_hostname\:$pubdir/SMV_Verification_Guide.pdf        $GUIDE_DIR/.
