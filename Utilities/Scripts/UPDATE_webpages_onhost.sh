#!/bin/bash

directory=$1
host=$2

ssh -q $host \( cd \~/$directory \; git checkout nist-pages \; git remote update \; git merge origin/nist-pages \)
ssh -q $host \( cd \~/$directory \; git describe --dirty --long \)
