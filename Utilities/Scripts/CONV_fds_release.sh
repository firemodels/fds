#!/bin/bash

bundleinfo=./bundle_setup
wikify=./wikify.py


echo obtaining FDS release notes wiki from the repository
svn export --quiet --force https://fds-smv.googlecode.com/svn/wiki/FDS_Release_Notes.wiki $bundleinfo/fds.wiki
echo Converting the wiki to html format
$wikify -r $bundleinfo/fds.wiki > $bundleinfo/release_notes.htm

echo obtaining FDS MPI Linux release notes wiki from the repository
svn export --quiet --force https://fds-smv.googlecode.com/svn/wiki/Running_FDS_MPI_on_Linux.wiki $bundleinfo/fds_linux.wiki
echo Converting the FDS MPI Linux wiki to html format
$wikify -r $bundleinfo/fds_linux.wiki > $bundleinfo/fds_linux.html

echo obtaining FDS MPI OX X release notes wiki from the repository
svn export --quiet --force https://fds-smv.googlecode.com/svn/wiki/Running_FDS_MPI_on_OSX.wiki $bundleinfo/fds_osx.wiki
echo Converting the FDS MPI OSX wiki to html format
$wikify -r $bundleinfo/fds_osx.wiki > $bundleinfo/fds_osx.html
