#!/bin/csh -f

set directory=$1
set    branch=$2
set      host=$3

echo
echo Updating the GIT repository: $directory, branch: $branch on host: $host to the latest revision
echo
ssh -q $host \( cd \~/$directory \; git checkout $branch \; git remote update \; git merge origin/$branch \)
ssh -q $host \( cd \~/$directory \; git describe --dirty  \)
