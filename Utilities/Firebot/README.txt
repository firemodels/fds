=====================================
= Firebot configuration information =
=====================================

# Firebot
# FDS automatIc veRification and validation tEst bot
# Kristopher Overholt
# 7/2/2012

=========
= About =
=========

Firebot is an automatic verification and validation test bot that is run at a regular interval (nightly).
More details on the Firebot build stages can be found in the FDS Configuration Management Plan.

=========================
= Firebot files/scripts =
=========================

# firebot_linux.sh

    This is the primary Firebot automated test script used on Linux with a queueing system (e.g., TORQUE).
    This script is invoked via crontab (details below).

# firebot_mac.sh

    This is the primary Firebot automated test script used on Mac OS X.
    This script is invoked via crontab (details below).

# firebot_website.py:

    This Python script generates an html file showing the status of recent Firebot test/build cycles.
    This script writes out an html file that can be 'cat' to /var/www/html/firebot/index.html.
    This script is invoked via crontab (details below).

# firebot_wiki_publish.sh:

    This bash shell script checks the Firebot build status wiki page (output by firebot_wiki.py) for changes.
    If changes to the build status wiki page are detected (compared to the SVN copy),
    Then Firebot commits the new status page to the SVN repository under fds.firebot@gmail.com.
    This requires a checkout of the wiki pages to exist in /home/firebot/fds_wiki via this command:
    svn checkout https://fds-smv.googlecode.com/svn/wiki/ /home/firebot/fds_wiki --username fds.firebot
    If you are doing this for the first time, be sure to manually send a test commit to the repo so that the username/password gets stored.
    This script is invoked via crontab (details below).

# firebot_wiki.py

    This Python script generates a Google Code wiki markup file showing the status of recent Firebot test/build cycles.
    This script writes out a plaintext .wiki file that can be 'cat' to /home/firebot/firebot/wiki_output/Firebot_Build_Status.wiki
    This script is invoked via crontab (details below).

# /usr/local/bin/run-one (from https://launchpad.net/ubuntu/+source/run-one)

    This wrapper script allows the execution of only one instance of a command and its args at a time.
    It does this by using a flock mechanism. This is useful to prepend to Firebot's crontab for firebot.sh so that,
    if Firebot is stalled for a long period of time, then multiple Firebot scripts do not run over each other.

===========
= Crontab =
===========

------------------------------------------------------------------------------------

#### The following information is in the Linux (Blaze) firebot user's crontab: ####

------------------------------------------------------------------------------------

PATH=/bin:/usr/bin:/usr/local/bin:/home/firebot/firebot:$PATH
MAILTO=""

#  ==========================
#  = Firebot status outputs =
#  ==========================

# Update Firebot status website every 10 minutes
*/10 * * * * firebot_website.py > /var/www/html/firebot/index.html

# Update Firebot build status wiki page every 10 minutes (XX:01, XX:11, and so on)
1,11,21,31,41,51 * * * * firebot_wiki.py > /home/firebot/firebot/wiki_output/Firebot_Build_Status.wiki

# Check every 10 minutes (XX:02, XX:12, and so on) for changes to Firebot status page
# If it has changed, upload Firebot build status wiki page to Google Code
2,12,22,32,42,52 * * * * firebot_wiki_publish.sh

#  ========================
#  = Firebot build script =
#  ========================

# Run svn update at 9:50 PM to get latest verison of Firebot
50 21 * * * cd ~/firebot ; svn up

# Run Firebot at 9:56 PM every night
# If no SVN argument is specified, then the latest SVN revision is used
# The run-once script maintains a lock to prevent the script from running twice
56 21 * * * run-one bash -lc firebot_linux.sh

------------------------------------------------------------------------------------

#### The following information is in the Mac (Bluesky) firebot user's crontab: ####

------------------------------------------------------------------------------------

PATH=/bin:/usr/bin:/usr/local/bin:/Users/firebot/firebot:$PATH
MAILTO=""

#  ========================
#  = Firebot build script =
#  ========================

# Run svn update at 9:50 PM to get latest verison of Firebot
50 21 * * * cd ~/firebot ; svn up

# Run Firebot at 9:56 PM every night
# If no SVN argument is specified, then the latest SVN revision is used
56 21 * * * cd ~/firebot ; ./firebot_mac.sh

------------------------------------------------------------------------------------
