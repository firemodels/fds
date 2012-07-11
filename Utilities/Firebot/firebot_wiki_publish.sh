#!/bin/bash

#  ========================
#  = Firebot wiki publish =
#  ========================
#
# Firebot
# FDS automatIc veRification and validation tEst bot
# Kristopher Overholt
# 7/2/2012
#
# This bash script can be run at a regular interval to commit
# the Firebot status wiki page to the SVN repository.

#  ===================
#  = Input variables =
#  ===================

SVNROOT_WIKI="/home/firebot/fds_wiki"
FIREBOT_DIR="/home/firebot/firebot"

#  ============================
#  = Primary script execution =
#  ============================

cd $FIREBOT_DIR

# Compare two wiki files, check for changes in Firebot status
   if [[ `diff ${FIREBOT_DIR}/wiki_output/Firebot_Build_Status.wiki ${SVNROOT_WIKI}/Firebot_Build_Status.wiki` == "" ]]
   then
      # No changes, exit script
      exit
   else
      cd $SVNROOT_WIKI
      svn update
      cp ${FIREBOT_DIR}/wiki_output/Firebot_Build_Status.wiki ${SVNROOT_WIKI}/Firebot_Build_Status.wiki
      svn commit -m 'Firebot: Update Firebot build status wiki page (Firebot_Build_Status)'
   fi
