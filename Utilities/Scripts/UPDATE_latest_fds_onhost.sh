#!/bin/bash

directory=$1
host=$2

ssh -q $host \( cd \~/$directory \; git checkout master \; git remote update \; git merge origin/master \; git merge firemodels/master  \)
ssh -q $host \( cd \~/$directory \; git describe --dirty --long \)
