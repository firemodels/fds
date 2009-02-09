#!/bin/csh -f 

# get test cases for Linux/OSX
# dummy change to force commit

set REPOS=~/FDS-SMV

# -------------- should not need to edit below ------------

set revision=$1

cd $REPOS/Utilities/Scripts/to_google/

set testdir=fds_test_cases_$revision

if (-d Test_cases) then
echo removing Test_cases
rm -rf Test_cases
endif
rm -rf $testdir.tar
rm -rf $testdir.tar.gz

svn export https://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Test_cases Test_cases
tar cvf $testdir.tar Test_cases/.
gzip $testdir.tar
