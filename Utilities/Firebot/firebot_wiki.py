#!/usr/bin/python

import os, re, time

#  =======================
#  = Firebot wiki status =
#  =======================
#
# Firebot
# FDS automatIc veRification and validation tEst bot
# Kristopher Overholt
# 7/2/2012
#
# This Python script can be 'cat' to ./wiki_output/firebot_status.wiki
# at a regular interval to display the status of the 100 most
# recent Firebot test runs on the Google Code Firebot status wiki page.

# Absolute path to Firebot history directory
firebot_history_dir = "/home/firebot/firebot/history"

# Print header information (doctype, head css, start of body)
print """#summary Displays FDS-SMV build status for recent build/test cycles.

=FDS Automatic Verification and Validation Test Bot=

----

==Firebot Build Status==

This page displays the build status for up to 100 of the most recent build/test cycles.

For additional information on Firebot, please refer to the FDS Configuration Management Plan.

===Build stages:===

  * Stage 1:  SVN operations
  * Stage 2a: Compile FDS DB
  * Stage 2b: Compile FDS MPI DB
  * Stage 3:  Run verification cases (short run)
  * Stage 4a: Compile FDS release
  * Stage 4b: Compile FDS MPI release
  * Stage 5:  Run verification cases (long run)
  * Stage 6a: Compile SMV utilities
  * Stage 6b: Compile SMV DB
  * Stage 6c: Make SMV pictures (debug mode)
  * Stage 6d: Compile SMV release
  * Stage 6e: Make SMV pictures (release mode)
  * Stage 6f: Make FDS pictures
  * Stage 7:  Matlab plotting and statistics
  * Stage 8:  Build FDS-SMV Guides
"""

print "----"

# Initialize list of history files
revision_files = os.listdir(firebot_history_dir)
revision_list = []

# Loop through history files and build list of unique revision numbers
for i in revision_files:
    rev_number_only = re.findall(r'\b\d+\b', i)
    if rev_number_only != []:
        rev_number_only = rev_number_only[0]
        revision_list.append(rev_number_only)

#  =======================================
#  = Print last 100 build status updates =
#  =======================================
for rev in sorted(revision_list, reverse=True)[:100]:

    #  ====================
    #  = Case: build fail =
    #  ====================
    if os.path.exists(firebot_history_dir + "/" + rev + "_errors.txt"):
        # Open error logs
        f = open (firebot_history_dir + "/" + rev + "_errors.txt","r")
        error_log = f.readlines()
        f.close()
        f_stage = open (firebot_history_dir + "/" + rev + ".txt","r")
        stage_log = f_stage.readlines()
        f_stage.close()

        # Write "Build fail" and hover of error log
        error_line = """Revision %(revision_num)s: [<font color="lightcoral"><b>x</b></font>] Build failure %(stage)s"""

        print
        print error_line % {'revision_num': rev, 'stage': (stage_log[0].split("failure")[1].strip())}

    #  ======================================
    #  = Case: build success, with warnings =
    #  ======================================
    elif os.path.exists(firebot_history_dir + "/" + rev + "_warnings.txt"):
        # Open warning log
        f = open (firebot_history_dir + "/" + rev + "_warnings.txt","r")
        warning_log = f.readlines()
        f.close()

        # Write "Build success, with warnings" and hover of warning log
        error_line = """Revision %(revision_num)s: [<font color="yellow"><b>-</b></font>] Build success, with warnings."""

        print
        print error_line % {'revision_num': rev}

    #  ========================
    #  = Case: build success! =
    #  ========================
    else:
        print
        print "Revision " + rev + ": [<font color='lime'><b>+</b></font>] Build success!"

    # Print time associated with history file
    print "<br>_<font color='#B0B0B0'>" + time.ctime(os.path.getmtime(firebot_history_dir + "/" + rev + ".txt")) + "</font>_"
