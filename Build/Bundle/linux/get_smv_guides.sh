#!/bin/bash

if [ -e ~/FDS_SMV_ENVpc.sh ]; then
  source ~/FDS_SMV_ENVpc.sh
else
  source ~/FDS_SMV_ENV.sh
fi

if [ ! -e ~/FDS_Guides ]; then
  mkdir ~/FDS_Guides
fi

pubdir=$firebotrepo/smv/Manuals
cp $pubdir/SMV_Technical_Reference_Guide/SMV_Technical_Reference_Guide.pdf ~/FDS_Guides/.
cp $pubdir/SMV_User__Guide/SMV_User_Guide.pdf                              ~/FDS_Guides/.
cp $pubdir/SMV_Verification_Guide/SMV_Verification_Guide.pdf               ~/FDS_Guides/.
