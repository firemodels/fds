#!/bin/bash
cat << EOF
<!DOCTYPE html>
<html><head><title>Firebot Build Status</title>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">

</head>
<body>
<h2>FDS Automatic Verification Summary Page</h2>

<h3>FDS Manuals</h3>
EOF

curdir=`pwd`
cd /var/www/html/firebot/manuals
ls -l *.pdf | awk '{printf("<li><a href=\"http://blaze.nist.gov/firebot/manuals/%s\">%s</a> <em>%s %s %s</em>\n",$9,$9,$6,$7,$8)}'
cd $curdir

cat << EOF
<h3>Firebot Status</h3>

This page displays the status for up to 30 of the most recent build/test cycles.<br>

EOF
