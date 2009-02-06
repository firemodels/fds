#!/bin/csh -f 

# get test cases for Linux/OSX

set REPOS=~/FDS-SMV

# -------------- should not need to edit below ------------

set revision=$1

cd $REPOS/Utilities/Scripts/to_google/

set testdir=FDS_Test_cases_$revision

if (-d $testdir) then
echo removing $testdir
rm -rf $testdir
endif
rm -rf $testdir.tar
rm -rf $testdir.tar.gz

svn export https://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Test_cases $testdir
cd $testdir
tar cvf ../$testdir.tar .
cd ..
gzip $testdir.tar
