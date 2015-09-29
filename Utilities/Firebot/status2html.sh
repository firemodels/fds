#!/bin/bash
CURDIR=`pwd`
cd history
listin=/tmp/list.in.$$
ls -tl *-????????.txt | head -30 |  awk '{system("head "  $9)}' | sort -t ';' -r -n -k 7 > $listin
cat $listin | awk -F ';' '{font="<font color=\"#00FF00\">";if($8=="2")font="<font color=\"#FF00FF\">";if($8=="3")font="<font color=\"#FF0000\">";printf("<p><a href=\"https://github.com/firemodels/fds-smv/commit/%s\">Revision: %s</a> %s%s</font><br>\n%s\n", $4,$5,font,$1,$2);}'
rm $listin
cd $CURDIR
