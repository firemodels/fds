#!/bin/bash
CURDIR=`pwd`
cpufrom=~/.smokebot/smv_times.csv

cat << EOF
<!DOCTYPE html>
<html><head><title>Smokebot Build Status</title>
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
          title: 'Smokebot CPU Time History',
          curveType: 'line',
          legend: { position: 'right' },
          colors: ['black'],
          pointSize: 5,
          hAxis:{ title: 'Day number'},
          vAxis:{ title: 'CPU Time (s)'}
        };
        options.legend = 'none';

        var chart = new google.visualization.LineChart(document.getElementById('curve_chart'));

        chart.draw(data, options);
      }
    </script>

</head>
<body>
<h2>Smokeview Automatic Verification Summary Page</h2>


<hr align='left'>

<div id="curve_chart" style="width: 500px; height: 300px"></div>
<h3>FDS/Smokeview Manuals</h3>
<a href="http://goo.gl/n1Q3WH">Manuals</a>

<h3>Smokebot Status</h3>

This page displays the status for up to 30 of the most recent build/test cycles.<br>
EOF

cd $CURDIR
