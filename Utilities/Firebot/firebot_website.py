#!/usr/bin/python

import os, re, time

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
# at a regular interval to display the status of the 100 most
# recent Firebot test runs.

# Absolute path to Firebot history directory
firebot_history_dir = "/home/firebot/firebot/history"

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

hr {
  width: 50%;
  height: 2px;
  color: #000;
  background-color: #000;
}
</style>

</head>
<body>
<h2><u>F</u><font color="#dddddd">DS automat</font><u>I</u><font color="#dddddd">c ve</font><u>R</u><font color="#dddddd">ification and validation t</font><u>E</u><font color="#dddddd">st</font> BOT</h2>

<hr align='left'>

<h3>Firebot Build Status</h3>

This page displays the build status for up to 100 of the most recent builds.<br>
To view a detailed error log, mouse over a fail or warning status.<br>

<br>
<span class="dropt">View a list of build stages<br>
<span style="width:700px;">
<b>Stage 1: </b> SVN operations <br><br>
<b>Stage 2a:</b> Compile FDS DB <br><br>
<b>Stage 2b:</b> Compile FDS MPI DB <br><br>
<b>Stage 3: </b> Run verification cases (short run) <br><br>
<b>Stage 4a:</b> Compile FDS release <br><br>
<b>Stage 4b:</b> Compile FDS MPI release <br><br>
<b>Stage 5: </b> Run verification cases (long run) <br><br>
<b>Stage 6a:</b> Compile SMV utilities <br><br>
<b>Stage 6b:</b> Compile SMV DB <br><br>
<b>Stage 6c:</b> Make SMV pictures (debug mode) <br><br>
<b>Stage 6d:</b> Compile SMV test <br><br>
<b>Stage 6e:</b> Make SMV pictures (test mode) <br><br>
<b>Stage 6f:</b> Compile SMV release <br><br>
<b>Stage 6g:</b> Make SMV pictures (release mode) <br><br>
<b>Stage 6h:</b> Make FDS pictures <br><br>
<b>Stage 7: </b> Matlab plotting and statistics <br><br>
<b>Stage 8: </b> Build FDS-SMV Guides <br>
</span></span>"""

current_time = time.ctime()
print "<br>"
print "<i>Status last updated: %s</i>" % (current_time)
print "<br><br><hr align='left'><br>"
print """Latest FDS Manuals: <a href="manuals/FDS_User_Guide.pdf">User Guide</a>, <a href="manuals/FDS_Technical_Reference_Guide.pdf">Technical Reference Guide</a>, <a href="manuals/FDS_Verification_Guide.pdf">Verification Guide</a>, <a href="manuals/FDS_Validation_Guide.pdf">Validation Guide</a><br><br>"""
print """Latest SMV Manuals: <a href="manuals/SMV_User_Guide.pdf">User Guide</a>, <a href="manuals/SMV_Verification_Guide.pdf">Verification Guide</a>"""
print "<br><br><hr align='left'><br>"

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
        error_line = """Revision number <b>%(revision_num)s</b>:
        <font style='BACKGROUND-COLOR: lightcoral'>[&#10007;]</font>
        <span class="dropt">Build failure<span style="width:700px;">"""

        print error_line % {'revision_num': rev}

        for line in error_log:
            print line + "<br>"

        print "</span>%s</span><br>" % (stage_log[0].split("failure")[1])

    #  ======================================
    #  = Case: build success, with warnings =
    #  ======================================
    elif os.path.exists(firebot_history_dir + "/" + rev + "_warnings.txt"):
        # Open warning log
        f = open (firebot_history_dir + "/" + rev + "_warnings.txt","r")
        warning_log = f.readlines()
        f.close()

        # Write "Build success, with warnings" and hover of warning log
        error_line = """Revision number <b>%(revision_num)s</b>:
        <font style='BACKGROUND-COLOR: yellow'>[&#8211;]</font>
        <span class="dropt">Build success, with warnings<br>
        <span style="width:700px;">"""
        print error_line % {'revision_num': rev}

        for line in warning_log:
            print line + "<br>"

        print "</span></span>"

    #  ========================
    #  = Case: build success! =
    #  ========================
    else:
        print "Revision number <b>" + rev + "</b>: <font style='BACKGROUND-COLOR: lime'>[&#10003;]</font> Build success!<br/>"

    # Print time associated with history file
    print "<font color='#B0B0B0'><i>" + time.ctime(os.path.getmtime(firebot_history_dir + "/" + rev + ".txt"))  + "</i></font>"
    print "<br><br>"

print """
</body>
</html>"""