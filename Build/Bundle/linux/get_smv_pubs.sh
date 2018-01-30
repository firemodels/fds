#!/bin/bash

source ../scripts/GET_ENV.sh

mkdir -p $GUIDE_DIR
pubdir=$smokebothome/.smokebot/pubs
cp $pubdir/SMV_Technical_Reference_Guide.pdf $GUIDE_DIR/.
cp $pubdir/SMV_User_Guide.pdf                $GUIDE_DIR/.
cp $pubdir/SMV_Verification_Guide.pdf        $GUIDE_DIR/.
