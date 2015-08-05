#!/bin/bash
headers=$1
CURDIR=`pwd`
cd history
listin=/tmp/list.in.$$
ls -tl *-????????.txt | head -30 |  awk '{system("head "  $9)}' > $listin
if [ "$headers" == "" ] ;then
cat << EOF
<!DOCTYPE html>
<html><head><title>Firebot Build Status</title>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">

</head>
<body>
<h2>FDS Automatic Verification Summary Page</h2>

<hr align='left'>

<h3>Firebot Status</h3>

This page displays the status for up to 30 of the most recent build/test cycles.<br>

EOF
fi
cat $listin | awk -F ';' '{printf "<p><a href=\"https://github.com/firemodels/fds-smv/commit/%s\">Revision: %s</a> %s<br>\n%s\n", $4,$5,$1,$2}'
if [ "$headers" == "" ] ;then
cat << EOF
<br><br>
<hr align='left'><br>
<i>Status last updated: `date`</i><br><br><br>

</body>
</html>
EOF
fi
rm $listin
cd $CURDIR
