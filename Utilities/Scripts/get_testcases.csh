#!/bin/csh -f 

# get test cases for Linux/OSX

set revision=$1
set REPOS=$2

cd $REPOS/Utilities/to_google/

set testdir=fds_examples_$revision

if (-d Examples) then
echo removing Examples
rm -rf Examples
endif
rm -rf $testdir.tar
rm -rf $testdir.tar.gz

svn export https://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Verification Examples
rm -rf Examples/Decaying_Isotropic_Turbulence
tar cvf $testdir.tar Examples/.
gzip $testdir.tar
