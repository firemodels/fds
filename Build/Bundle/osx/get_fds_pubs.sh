#!/bin/bash

source ../scripts/GET_ENV.sh

mkdir -p $GUIDE_DIR
pubdir=$firebothome/.firebot/pubs
scp -q $linux_hostname\:$pubdir/FDS_Config_Management_Plan.pdf    $GUIDE_DIR/.
scp -q $linux_hostname\:$pubdir/FDS_Technical_Reference_Guide.pdf $GUIDE_DIR/.
scp -q $linux_hostname\:$pubdir/FDS_User_Guide.pdf                $GUIDE_DIR/.
scp -q $linux_hostname\:$pubdir/FDS_Validation_Guide.pdf          $GUIDE_DIR/.
scp -q $linux_hostname\:$pubdir/FDS_Verification_Guide.pdf        $GUIDE_DIR/.
