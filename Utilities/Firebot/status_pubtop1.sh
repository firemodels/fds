#!/bin/bash
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
