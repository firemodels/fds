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
span.dropt:hover span {margin: 5px 0 0 170px; background: #ffffff; z-index:6;}

hr {
  width: 50%;
  height: 2px;
  color: #000;
  background-color: #000;
}
</style>

</head>
<body>
<h1>F<font size=2>DS automat</font>I<font size=2>c ve</font>R<font size=2>ification and validation t</font>E<font size=2>st</font> BOT</h1>

<hr align='left'>

<h3>Firebot Build Status</h3>

This page displays the build status for up to 100 of the most recent builds.</br>
Hover over a fail or warning status for an error log.</br>"""

current_time = time.ctime()
print "</br>"
print "<i>Status last updated: %s</i>" % (current_time)
print "<hr align='left'></br>"

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
        # Open error log
        f = open (firebot_history_dir + "/" + rev + "_errors.txt","r")
        error_log = f.readlines()
        f.close()

        # Write "Build fail" and hover of error log
        error_line = """Revision number <b>%(revision_num)s</b>:
        <font style='BACKGROUND-COLOR: lightcoral'>[&#10007;]</font>
        <span class="dropt">Build fail</br>
        <span style="width:700px;">"""
        print error_line % {'revision_num': rev}

        for line in error_log:
            print line + "</br>"

        print "</span></span>"

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
        <span class="dropt">Build success, with warnings</br>
        <span style="width:700px;">"""
        print error_line % {'revision_num': rev}

        for line in warning_log:
            print line + "</br>"

        print "</span></span>"

    #  ========================
    #  = Case: build success! =
    #  ========================
    else:
        print "Revision number <b>" + rev + "</b>: <font style='BACKGROUND-COLOR: lime'>[&#10003;]</font> Build success!<br/>"

    # Print time associated with history file
    print time.ctime(os.path.getmtime(firebot_history_dir + "/" + rev + ".txt"))
    print "</br></br>"

print """
</body>
</html>"""