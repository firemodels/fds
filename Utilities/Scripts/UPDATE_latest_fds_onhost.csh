#!/bin/csh -f

set directory=$1
set      host=$2

echo Updating the GIT repository $directory on $host to the latest revision
echo
ssh -q $host \( cd \~/$directory \; git checkout development \; git remote update \; git merge origin/development \; git merge firemodels/development  \)
