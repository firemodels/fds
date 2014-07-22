#!/bin/bash

bundleinfo=./bundle_setup
wikify=./wikify.py


echo obtaining FDS release notes wiki from the repository
svn export --quiet --force https://fds-smv.googlecode.com/svn/wiki/FDS_Release_Notes.wiki $bundleinfo/fds.wiki
echo Converting the wiki to html format
$wikify -r $bundleinfo/fds.wiki > $bundleinfo/release_notes.htm
