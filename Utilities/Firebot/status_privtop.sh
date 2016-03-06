#!/bin/bash
curdir=`pwd`
cpufrom=~/.firebot/fds_times.csv

cat << EOF
<!DOCTYPE html>
<html><head><title>Firebot Build Status</title>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
    <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
    <script type="text/javascript">
      google.charts.load('current', {'packages':['corechart']});
      google.charts.setOnLoadCallback(drawChart);

      function drawChart() {
        var data = google.visualization.arrayToDataTable([
          ['Days since Jan 1, 2016', 'CPU Time (s)'],
EOF

sort -n -k 1 -t , $cpufrom | tail -30 | awk -F ',' '{ printf("[%s,%s],\n",$1,$2) }'

cat << EOF
        ]);

        var options = {
          title: 'Firebot CPU History',
          curveType: 'line',
          legend: { position: 'bottom' },
          colors: ['black'],
          pointSize: 5
        };

        var chart = new google.visualization.LineChart(document.getElementById('curve_chart'));

        chart.draw(data, options);
      }
    </script>

</head>
<body>
<h2>FDS Automatic Verification Summary Page</h2>
<div id="curve_chart" style="width: 500px; height: 300px"></div>
EOF

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
<h3>Firebot Status</h3>

Firebot status for up to 30 previous build/test cycles.<br>
EOF
