=====================================
= Firebot configuration information =
=====================================

# Firebot
# FDS automatIc veRification and validation tEst bot
# Kristopher Overholt
# 7/2/2012

RUNNING FIREBOT: To run Firebot manually on Blaze at NIST, use the run_firebot_linux.sh script.

=========
= About =
=========

Firebot is an automatic verification and validation test bot that is run at a regular interval (nightly).
More details on the Firebot build stages can be found in the FDS Configuration Management Plan.

======================
= Installing Firebot =
======================

1. Create an account for Firebot on your machine. This account should be named "firebot" if possible.
   The use of a separate account is recommended so that Firebot has a clean and undisturbed copy of
   the repositories to work with.

2. In Firebot's newly created home directory, perform an SVN checkout on the Firebot portion of the
   FDS-SMV repository using:

      svn checkout https://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Utilities/Firebot firebot

3. cd to the newly created ~/firebot directory

4. Run the ./firebot_linux.sh or ./firebot_mac.sh command, then the automated Firebot build process
   will begin and will create directories and check out repositories as needed.

(Optional) When running Firebot automatically using the cron system, the use of the run-one script is recommended,
which is available from https://launchpad.net/ubuntu/+source/run-one. This script ensures that one instance of
Firebot will not run over another instance of Firebot, which would cause file conflicts. The run-one script
can be installed and used locally, and you can see examples of its usage in the crontab information below.

=========================
= Firebot files/scripts =
=========================

# run_firebot_linux.sh

    This script can be used to start Firebot manually. First, the script performs an "svn update".
    Then, the script runs Firebot using the run-one script and backgrounds Firebot.

# firebot_linux.sh

    This is the primary Firebot automated test script used on Linux with a queueing system (e.g., TORQUE).
    This script is invoked via crontab (details below).
    By default, Firebot is run in verification mode, but the -v option can be used to run Firebot in validation mode.

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

#### The following information is in the Linux firebot user's crontab: ####

------------------------------------------------------------------------------------

PATH=/bin:/usr/bin:/usr/local/bin:/home2/smokevis2/firebot/firebot:$PATH
MAILTO=""

#  ==========================
#  = Firebot status outputs =
#  ==========================

# Update Firebot status website every 10 minutes
*/10 * * * * firebot_website.py > /var/www/html/firebot/index.html

# Update Firebot build status wiki page every 10 minutes (XX:01, XX:11, and so on)
1,11,21,31,41,51 * * * * firebot_wiki.py > ~/firebot/wiki_output/Firebot_Build_Status.wiki

# Check every 10 minutes (XX:05, XX:15, and so on) for changes to Firebot status page
# If it has changed, upload Firebot build status wiki page to Google Code
5,15,25,35,45,55 * * * * firebot_wiki_publish.sh

#  ========================
#  = Firebot build script =
#  ========================

# Run svn update at 9:50 PM to get latest verison of Firebot
50 21 * * * cd ~/firebot ; svn revert * ; svn up

# Run Firebot at 9:56 PM every night
# If no SVN argument is specified, then the latest SVN revision is used
# The run-once script maintains a lock to prevent the script from running twice
56 21 * * * run-one bash -lc firebot_linux.sh

# ============================
# = DiskHog disk space alert =
# ============================

# Run DiskHog script daily at 8 AM to check disk usage on all nodes
# Sends email alert if any disk on any node is more than 90% full
00 08 * * * diskhog.sh

# =====================================
# = Torque job processing error alert =
# =====================================

# Run torque error checking script once a day at
# 23:55 to check the torque log for processing errors.
# Sends email alert if any "post job processing errors"
# are found in the daily torque log.
55 23 * * * chktorque.sh

------------------------------------------------------------------------------------

#### The following information is in the Linux validationbot user's crontab: ####

------------------------------------------------------------------------------------

PATH=/bin:/usr/bin:/usr/local/bin:/home2/smokevis2/validationbot/firebot:$PATH
MAILTO=""

#  ========================
#  = Validationbot script =
#  ========================

# Run svn update every hour at XX:25 to get latest verison of Firebot
25 * * * * cd ~/firebot ; svn revert * ; svn up

# Run Validationbot every hour at xx:30
# The run-once script maintains a lock to prevent the script from running twice
30 * * * * run-one bash -lc "firebot_linux.sh -s -v -y"

------------------------------------------------------------------------------------

#### The following information is in the Mac firebot user's crontab: ####

------------------------------------------------------------------------------------

PATH=/bin:/usr/bin:/usr/local/bin:/usr/sbin:/Users/firebot/firebot:$PATH
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
