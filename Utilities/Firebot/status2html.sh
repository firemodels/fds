#!/bin/bash
CURDIR=`pwd`
cd history
listin=/tmp/list.in.$$
ls -tl *-????????.txt | head -30 |  awk '{system("head "  $9)}' | sort -t ';' -r -n -k 7 > $listin
cat $listin | awk -F ';' '{font="<font color=\"#2EB82E\">";if($8=="2")font="<font color=\"#CCCC00\">";if($8=="3")font="<font color=\"#E61616\">";printf("<p><a href=\"https://github.com/firemodels/fds-smv/commit/%s\">Revision: %s</a> %s%s</font><br>\n%s\n", $4,$5,font,$1,$2);}'
rm $listin
cd $CURDIR
