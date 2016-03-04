#!/bin/bash
cat << EOF
<!DOCTYPE html>
<html><head><title>Firebot Build Status</title>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">

</head>
<body>
<h2>FDS Automatic Verification Summary Page</h2>
EOF

curdir=`pwd`
cd /var/www/html/firebot/manuals
cat << EOF
<h3>FDS Manuals</h3>
<ul>
EOF
ls -l *.pdf | awk '{printf("<li><a href=\"http://blaze.nist.gov/firebot/manuals/%s\">%s</a> <em>%s %s %s</em>\n",$9,$9,$6,$7,$8)}'
cat << EOF
</ul>
EOF
cd $curdir

cd /var/www/html/smokebot/manuals
cat << EOF
<h3>Smokeview Manuals</h3>
<ul>
EOF
ls -l *.pdf | grep -v geom|awk '{printf("<li><a href=\"http://blaze.nist.gov/smokebot/manuals/%s\">%s</a> <em>%s %s %s</em>\n",$9,$9,$6,$7,$8)}'
cat << EOF
</ul>
EOF
cd $curdir

cat << EOF
<h3>Benchmark Time History</h3>
<IMG width=500 SRC="fds_times.png"><br>
<h3>Firebot Status</h3>

Firebot status for up to 30 previous build/test cycles.<br>
EOF
