#!/bin/bash
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

<hr align='left'>

<div id="curve_chart" style="width: 500px; height: 300px"></div>
<h3>FDS/Smokeview Manuals</h3>
<a href="http://goo.gl/n1Q3WH">Manuals</a>

<h3>Firebot Status</h3>

This page displays the status for up to 30 of the most recent build/test cycles.<br>
EOF
