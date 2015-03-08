#!/usr/bin/python

import os
import re
import time

#  ==========================
#  = Firebot website status =
#  ==========================
#
# Firebot
# FDS automatIc veRification and validation tEst bot
# Kristopher Overholt
# 7/2/2012
#
# This Python script can be 'cat' to /var/www/html/firebot/index.html
# at a regular interval to display the status of the 30 most
# recent Firebot test runs.

# Absolute path to Firebot history directory
firebot_history_dir = "/home2/smokevis2/firebot/firebot/history"

# Print header information (doctype, head css, start of body)
print """<!DOCTYPE html>
<html><head><title>Firebot Build Status</title>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">

<style type="text/css">
span.dropt {border-bottom: thin dotted; background: #ffeedd;}
span.dropt:hover {text-decoration: none; background: #ffffff; z-index: 6; }
span.dropt span {position: absolute; left: -9999px;
  margin: 20px 0 0 0px; padding: 3px 3px 3px 3px;
  border-style:solid; border-color:black; border-width:1px; z-index: 6;}
span.dropt:hover span {left: 2%; background: #ffffff;}
span.dropt span {position: absolute; left: -9999px;
  margin: 4px 0 0 0px; padding: 3px 3px 3px 3px;
  border-style:solid; border-color:black; border-width:1px;}
span.dropt:hover span {margin: 5px 0 0 170px; background: #ebf7f8; z-index:6;}

/* Keep link color as unvisited */
a {
  color:#0000FF;
}

/* Set spacing between list items */
li {
   line-height: 1.5em;
}

/* Set primary font */
body{
margin-left: 50px;
font: 16px Arial, sans-serif;
}

hr {
  width: 65%;
  height: 2px;
  color: #000;
  background-color: #000;
}
</style>

</head>
<body>
<h2>FDS Automatic Verification and Validation Test Bot</h2>

<hr align='left'>

<h3>Firebot Build Status</h3>

This page displays the build status for up to 30 of the most recent build/test cycles.<br>
To view a detailed error log, mouse over a fail or warning status.<br>

<br>
<span class="dropt">View a list of build stages<br>
<span style="width:700px;">
<b>Stage 1:</b> SVN operations <br>
<b>Stage 2a:</b> Compile and inspect FDS debug <br>
<b>Stage 2b:</b> Compile FDS MPI debug <br>
<b>Stage 3:</b> Run verification cases (debug mode) <br>
<b>Stage 4a:</b> Compile FDS release <br>
<b>Stage 4b:</b> Compile FDS MPI release <br>
<b>Stage 5pre:</b> Compile SMV utilities <br>
<b>Stage 5:</b> Run verification cases (release mode) <br>
<b>Stage 6a:</b> Compile SMV debug <br>
<b>Stage 6b:</b> Make SMV pictures (debug mode) <br>
<b>Stage 6c:</b> Compile SMV release <br>
<b>Stage 6d:</b> Make SMV pictures (release mode) <br>
<b>Stage 6e:</b> Make FDS pictures <br>
<b>Stage 7a:</b> Matlab plotting and statistics (verification) <br>
<b>Stage 7b:</b> Matlab plotting (validation) <br>
<b>Stage 7c:</b> FDS run time statistics <br>
<b>Stage 8:</b> Build FDS-SMV guides <br>
</span></span>"""

print "<br><hr align='left'><br>"

print "FDS Manuals: <i>(last updated " + time.ctime(os.path.getmtime("/var/www/html/firebot/manuals/FDS_User_Guide.pdf")) + ")</i><br><br>"

print """<li><a href="manuals/FDS_User_Guide.pdf">FDS User Guide</a></li>
         <li><a href="manuals/FDS_Technical_Reference_Guide.pdf">FDS Technical Reference Guide</a></li>
         <li><a href="manuals/FDS_Verification_Guide.pdf">FDS Verification Guide</a></li>
         <li><a href="manuals/FDS_Validation_Guide.pdf">FDS Validation Guide</a></li>
         <li><a href="manuals/FDS_Configuration_Management_Plan.pdf">FDS Configuration Management Plan</a></li>
         <br>
         <li><a href="correlation_guide/Correlation_Guide.pdf">Correlation Guide</a></li>"""

print "<p>Smokeview Manuals: <i>(last updated " + time.ctime(os.path.getmtime("/var/www/html/smokebot/manuals/SMV_User_Guide.pdf")) + ")</i><br><br>"

print """<li><a href="../smokebot/manuals/SMV_User_Guide.pdf">SMV User Guide</a></li>
         <li><a href="../smokebot/manuals/SMV_Technical_Reference_Guide.pdf">SMV Technical Reference Guide</a></li>
         <li><a href="../smokebot/manuals/SMV_Verification_Guide.pdf">SMV Verification Guide</a></li>"""
print "<br><hr align='left'><br>"

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
#  = Print last 30 build status updates =
#  =======================================
for rev in sorted(revision_list, reverse=True)[:30]:

    #  ====================
    #  = Case: build fail =
    #  ====================
    if os.path.exists(firebot_history_dir + "/" + rev + "_errors.txt"):
        # Open error logs
        f = open(firebot_history_dir + "/" + rev + "_errors.txt", "r")
        error_log = f.readlines()
        f.close()
        f_stage = open(firebot_history_dir + "/" + rev + ".txt", "r")
        stage_log = f_stage.readlines()
        f_stage.close()

        # Write "Build fail" and hover of error log
        error_line = """<a href="http://code.google.com/p/fds-smv/source/detail?r=%(revision_num)s">Revision %(revision_num)s</a>:
        <font style='BACKGROUND-COLOR: lightcoral'>[&#10007;]</font>
        <span class="dropt">Build failure.<span style="width:700px;">"""

        print error_line % {'revision_num': rev}

        for line in error_log:
            print line + "<br>"

        print "</span></span><br>"

    #  ======================================
    #  = Case: build success, with warnings =
    #  ======================================
    elif os.path.exists(firebot_history_dir + "/" + rev + "_warnings.txt"):
        # Open warning log
        f = open(firebot_history_dir + "/" + rev + "_warnings.txt", "r")
        warning_log = f.readlines()
        f.close()

        # Write "Build success, with warnings" and hover of warning log
        error_line = """<a href="http://code.google.com/p/fds-smv/source/detail?r=%(revision_num)s">Revision %(revision_num)s</a>:
        <font style='BACKGROUND-COLOR: yellow'>[&#8211;]</font>
        <span class="dropt">Build success, with warnings.<br>
        <span style="width:700px;">"""
        print error_line % {'revision_num': rev}

        for line in warning_log:
            print line + "<br>"

        print "</span></span>"

    #  ========================
    #  = Case: build success! =
    #  ========================
    else:
        print "<a href='http://code.google.com/p/fds-smv/source/detail?r=" + rev + "'>Revision " + rev + "</a>: <font style='BACKGROUND-COLOR: lime'>[&#10003;]</font> Build success!<br/>"

    # Print time associated with history file
    print "<font color='#B0B0B0'><i>" + time.ctime(os.path.getmtime(firebot_history_dir + "/" + rev + ".txt")) + "</i></font>"
    print "<br><br>"

print "<hr align='left'><br>"

current_time = time.ctime()
print "<i>Status last updated %s</i><br><br><br>" % (current_time)

print """
</body>
</html>"""
