#!/bin/bash

bundleinfo=./bundle_setup
wikify=./wikify.py


echo Converting the FDS release notes from wiki to html format
$wikify -r $bundleinfo/FDS_Release_Notes.wiki > $bundleinfo/release_notes.htm
