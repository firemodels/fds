#!/bin/bash
historydir=~/.firebot/history
BODY=
TITLE=Firebot
SOPT=

while getopts 'bs' OPTION
do
case $OPTION  in
  b)
   BODY="1"
   ;;
  s)
   historydir=~/.smokebot/history
   TITLE=Smokebot
   SOPT=-s
   ;;
esac
done
shift $(($OPTIND-1))


if [ "$BODY" == "" ]; then
cat << EOF
<!DOCTYPE html>
<html><head><title>$TITLE Build Status</title>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
    <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
    <script type="text/javascript">
      google.charts.load('current', {'packages':['corechart']});
      google.charts.setOnLoadCallback(drawChart);

      function drawChart() {
        var data = google.visualization.arrayToDataTable([
          ['Days since Jan 1, 2016', 'Benchmark Time (s)'],
EOF

./make_timelist.sh $SOPT | sort -n -k 1 -t , | tail -30 | awk -F ',' '{ printf("[%s,%s],\n",$1,$2) }'

cat << EOF
        ]);

        var options = {
          title: '$TITLE Time History',
          curveType: 'line',
          legend: { position: 'right' },
          colors: ['black'],
          pointSize: 5,
          hAxis:{ title: 'Day'},
          vAxis:{ title: 'Benchmark Time (s)'}
        };
        options.legend = 'none';

        var chart = new google.visualization.LineChart(document.getElementById('curve_chart'));

        chart.draw(data, options);
      }
    </script>

</head>
<body>
<h2>$TITLE Summary Page</h2>

<hr align='left'>

<div id="curve_chart" style="width: 500px; height: 300px"></div>
<h3>FDS/Smokeview Manuals</h3>
<a href="http://goo.gl/n1Q3WH">Manuals</a>

<h3>Status</h3>

EOF
fi

CURDIR=`pwd`
listin=/tmp/list.in.$$
cd $historydir
ls -tl *-????????.txt | awk '{system("head "  $9)}' | sort -t ';' -r -n -k 7 > $listin
cat $listin | head -30 | \
             awk -F ';' '{cputime="Benchmark time: "$9" s";\
                          if($9=="")cputime="";\
                          font="<font color=\"#00FF00\">";\
                          if($8=="2")font="<font color=\"#FF00FF\">";\
                          if($8=="3")font="<font color=\"#FF0000\">";\
                          printf("<p><a href=\"https://github.com/firemodels/fds-smv/commit/%s\">Revision: %s</a>%s %s</font><br>\n",$4,$5,font,$1);\
                          printf("Revision date: %s\n",$2);\ 
                          if($9!="")printf("%s <br>\n",cputime);}'
rm $listin

if [ "$BODY" == "" ]; then
cat << EOF
<br><br>
<hr align='left'><br>
<i>Updated: `date`</i><br><br><br>

</body>
</html>
EOF
fi

cd $CURDIR
