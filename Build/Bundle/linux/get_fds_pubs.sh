#!/bin/bash

if [ -e ~/FDS_SMV_ENVpc.sh ]; then
  source ~/FDS_SMV_ENVpc.sh
else
  source ~/FDS_SMV_ENV.sh
fi

if [ ! -e ~/FDS_Guides ]; then
  mkdir ~/FDS_Guides
fi

pubdir=$firebothome/.firebot/pubs
cp $pubdir/FDS_Config_Management_Plan.pdf    ~/FDS_Guides/.
cp $pubdir/FDS_Technical_Reference_Guide.pdf ~/FDS_Guides/.
cp $pubdir/FDS_User_Guide.pdf                ~/FDS_Guides/.
cp $pubdir/FDS_Validation_Guide.pdf          ~/FDS_Guides/.
cp $pubdir/FDS_Verification_Guide.pdf        ~/FDS_Guides/.
