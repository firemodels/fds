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

SVNROOT_WIKI="/home2/smokevis2/firebot/fds_wiki"
FIREBOT_DIR="/home2/smokevis2/firebot/firebot"

#  ============================
#  = Primary script execution =
#  ============================

cd $SVNROOT_WIKI
svn cleanup
svn update

# Compare two wiki files, check for changes in Firebot status
   if [[ `diff ${FIREBOT_DIR}/wiki_output/Firebot_Build_Status.wiki ${SVNROOT_WIKI}/Firebot_Build_Status.wiki` == "" ]]
   then
      # No changes, exit script
      exit
   else
      # Changes detected, copy and commit new version
      cp ${FIREBOT_DIR}/wiki_output/Firebot_Build_Status.wiki ${SVNROOT_WIKI}/Firebot_Build_Status.wiki
      svn commit -m 'Firebot: Update Firebot build status wiki page (Firebot_Build_Status)'
   fi
