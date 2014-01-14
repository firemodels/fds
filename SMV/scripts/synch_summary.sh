#!/bin/bash
FROMDIR=/home2/firevis2/smokebot/FDS-SMV/Manuals/SMV_Summary
TODIR=~/gfiles/SMV_Summary
touch ~/bin/synch_start
rsync -avzu --exclude=.svn/ --delete $FROMDIR/ $TODIR
cd ~/gfiles
grive
touch ~/bin/synch_end
