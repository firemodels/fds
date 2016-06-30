#!/bin/bash
CURDIR=`pwd`

cat << EOF
<html>
<head>
    <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
    <script type="text/javascript">
      google.charts.load('current', {'packages':['corechart']});
      google.charts.setOnLoadCallback(drawChart);

      function drawChart() {
        var data = google.visualization.arrayToDataTable([
          ['Days since Jan 1, 2016', 'Benchmark Time (s)'],
EOF

./make_timelist.sh -s | sort -n -k 1 -t , | tail -30 | awk -F ',' '{ printf("[%s,%s],\n",$1,$2) }'

cat << EOF
        ]);

        var options = {
          title: 'Smokebot Time History',
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
EOF
cat << EOF
<TITLE>Smokeview Verification Test</TITLE>
</HEAD>
<BODY BGCOLOR="#FFFFFF" >
<h2>Smokeview Verification Tests</h2>
<p>
<a href="movies.html"><font size=5>Animations</font></a> -
<font size=5>Stills</font> -  
<a href="manuals.html"><font size=5>Manuals</font></a> 

<p><table>
<tr><td valign=middle>Version info:</td><td><img width=400 src="images/fds_smv_version.png"></td></tr>
</table>
<div id="curve_chart" style="width: 500px; height: 300px"></div>

<ul>
<li><a href="#slice">Slice files</a>
<li><a href="#boundary">Boundary files</a>
<li><a href="#isosurfaces">Isosurfaces</a>
<li><a href="#particles">Particles</a>
<li><a href="#3dsmoke">3D Smoke</a>
<li><a href="#plot3d">PLOT3D</a>
<li><a href="#zone">Zone fire</a>
<li><a href="#wui">WUI</a>
<li><a href="#obstacles">Obstacles</a>
<li><a href="#vents">Vents</a>
<li><a href="#stereo">Stereo</a>
<li><a href="#misc">Miscellaneous</a>
<li><a href="#obj2state">Objects - two states</a>
<li><a href="#obj1state">Objects - one state</a>
<li><a href="#objelem">Objects - elementary</a>
</ul>

<p><hr>

<a name="slice"></a>
<p><table>
<tr><th></th><th colspan=3> Slice</th></tr>
<tr>
<tr>
<th>node centered</th>
<td><a href="images/plume5c_slice_00.png"><img width=200 src="images/plume5c_slice_00.png"></a></td>
<td><a href="images/plume5c_slice_10.png"><img width=200 src="images/plume5c_slice_10.png"></a></td>
<td><a href="images/plume5c_slice_30.png"><img width=200 src="images/plume5c_slice_30.png"></a></td>
</tr>

<tr>
<th>cell centered</th>
<td><a href="images/plume5c_slice_cell_00.png"><img width=200 src="images/plume5c_slice_cell_00.png"></a></td>
<td><a href="images/plume5c_slice_cell_10.png"><img width=200 src="images/plume5c_slice_cell_10.png"></a></td>
<td><a href="images/plume5c_slice_cell_30.png"><img width=200 src="images/plume5c_slice_cell_30.png"></a></td>
</tr>

<tr>
<th>node centered<br>chopped</th>
<td><a href="images/plume5c_slice_chop_00.png"><img width=200 src="images/plume5c_slice_chop_00.png"></a></td>
<td><a href="images/plume5c_slice_chop_10.png"><img width=200 src="images/plume5c_slice_chop_10.png"></a></td>
<td><a href="images/plume5c_slice_chop_30.png"><img width=200 src="images/plume5c_slice_chop_30.png"></a></td>
</tr>

<tr>
<th>cell centered<br>chopped</th>
<td><a href="images/plume5c_slice_cellchop_00.png"><img width=200 src="images/plume5c_slice_cellchop_00.png"></a></td>
<td><a href="images/plume5c_slice_cellchop_10.png"><img width=200 src="images/plume5c_slice_cellchop_10.png"></a></td>
<td><a href="images/plume5c_slice_cellchop_30.png"><img width=200 src="images/plume5c_slice_cellchop_30.png"></a></td>
</tr>

<tr>
<th>vector</th>
<td><a href="images/plume5c_vslice_00.png"><img width=200 src="images/plume5c_vslice_00.png"></a></td>
<td><a href="images/plume5c_vslice_10.png"><img width=200 src="images/plume5c_vslice_10.png"></a></td>
<td><a href="images/plume5c_vslice_30.png"><img width=200 src="images/plume5c_vslice_30.png"></a></td>
</tr>

<tr>
<th>general orientation</th>
<td><a href="images/plume5c_gslice_00.png"><img width=200 src="images/plume5c_gslice_00.png"></a></td>
<td><a href="images/plume5c_gslice_10.png"><img width=200 src="images/plume5c_gslice_10.png"></a></td>
<td><a href="images/plume5c_gslice_30.png"><img width=200 src="images/plume5c_gslice_30.png"></a></td>
</tr>

<tr>
<th>differenced</th>
<td><a href="images/plume5c_diff_slice_00.png"><img width=200 src="images/plume5c_diff_slice_00.png"></a></td>
<td><a href="images/plume5c_diff_slice_10.png"><img width=200 src="images/plume5c_diff_slice_10.png"></a></td>
<td><a href="images/plume5c_diff_slice_30.png"><img width=200 src="images/plume5c_diff_slice_30.png"></a></td>
</tr>

<tr>
<th>coordinates</th>
<td><a href="images/cell_test_3D.png"><img width=200 src="images/cell_test_3D.png"></a></td>
<td><a href="images/cell_test_x.png"><img width=200 src="images/cell_test_x.png"></a></td>
<td><a href="images/cell_test_y.png"><img width=200 src="images/cell_test_y.png"></a></td>
<td><a href="images/cell_test_z.png"><img width=200 src="images/cell_test_z.png"></a></td>
</tr>

<tr>
<th>FED</th>
<td><a href="images/fed_test_030.png"><img width=200 src="images/fed_test_030.png"></a></td>
<td><a href="images/fed_test_060.png"><img width=200 src="images/fed_test_060.png"></a></td>
<td><a href="images/fed_test_120.png"><img width=200 src="images/fed_test_120.png"></a></td>
<td><a href="images/fed_test_360.png"><img width=200 src="images/fed_test_360.png"></a></td>
</tr>

<tr align=center>
<th>cell mask</th>
<td><a href="images/slicemask_cell1.png"><img width=200 src="images/slicemask_cell1.png"></a></td>
<td><a href="images/slicemask_cell2.png"><img width=200 src="images/slicemask_cell2.png"></a></td>
<td><a href="images/slicemask_cell3.png"><img width=200 src="images/slicemask_cell3.png"></a></td>
</tr>
<tr align=center><td></td><td>below</td><td>through</td><td>above</td></tr>

<tr align=center>
<th>node mask</th>
<td><a href="images/slicemask_node1.png"><img width=200 src="images/slicemask_node1.png"></a></td>
<td><a href="images/slicemask_node2.png"><img width=200 src="images/slicemask_node2.png"></a></td>
<td><a href="images/slicemask_node3.png"><img width=200 src="images/slicemask_node3.png"></a></td>
<td><a href="images/slicemask_node4.png"><img width=200 src="images/slicemask_node4.png"></a></td>
</tr>
<tr align=center><td></td><td>below</td><td>below on boundary</td><td>above on boundary</td><td>above</td></tr>
</table>

<p><hr>


<a name="boundary"></a>
<p><table>
<tr><th></th><th colspan=3>Boundary</th></tr>
<tr>
<th>node centered</th>
<td><a href="images/plume5c_bound_00.png"><img width=200 src="images/plume5c_bound_00.png"></a></td>
<td><a href="images/plume5c_bound_10.png"><img width=200 src="images/plume5c_bound_10.png"></a></td>
<td><a href="images/plume5c_bound_30.png"><img width=200 src="images/plume5c_bound_30.png"></a></td>
</tr>

<tr>
<th>cell centered</th>
<td><a href="images/plume5c_bound_cell_00.png"><img width=200 src="images/plume5c_bound_cell_00.png"></a></td>
<td><a href="images/plume5c_bound_cell_10.png"><img width=200 src="images/plume5c_bound_cell_10.png"></a></td>
<td><a href="images/plume5c_bound_cell_30.png"><img width=200 src="images/plume5c_bound_cell_30.png"></a></td>
</tr>
</table>

<p><hr>
<a name="isosurfaces"></a>
<p><table>
<tr><th></th><th colspan=3> Isosurfaces</th></tr>
<tr>
<th>points</th>
<td><a href="images/plume5c_part_01.png"><img width=200 src="images/plume5c_iso_points_00.png"></a></td>
<td><a href="images/plume5c_part_10.png"><img width=200 src="images/plume5c_iso_points_10.png"></a></td>
<td><a href="images/plume5c_part_30.png"><img width=200 src="images/plume5c_iso_points_30.png"></a></td>
</tr>

<tr>
<th>outline</th>
<td><a href="images/plume5c_part_01.png"><img width=200 src="images/plume5c_iso_outline_00.png"></a></td>
<td><a href="images/plume5c_part_10.png"><img width=200 src="images/plume5c_iso_outline_10.png"></a></td>
<td><a href="images/plume5c_part_30.png"><img width=200 src="images/plume5c_iso_outline_30.png"></a></td>
</tr>

<tr>
<th>solid</th>
<td><a href="images/plume5c_part_01.png"><img width=200 src="images/plume5c_iso_solid_00.png"></a></td>
<td><a href="images/plume5c_part_10.png"><img width=200 src="images/plume5c_iso_solid_10.png"></a></td>
<td><a href="images/plume5c_part_30.png"><img width=200 src="images/plume5c_iso_solid_30.png"></a></td>
</tr>
</table>

<p><hr>
<a name="particles"></a>
<p><table>
<tr><th></th><th colspan=3> Particles</th></tr>
<tr>
<th>points</th>
<td><a href="images/plume5c_part_01.png"><img width=200 src="images/plume5c_part_01.png"></a></td>
<td><a href="images/plume5c_part_10.png"><img width=200 src="images/plume5c_part_10.png"></a></td>
<td><a href="images/plume5c_part_30.png"><img width=200 src="images/plume5c_part_30.png"></a></td>
</tr>

<tr>
<th>streaks</th>
<td><a href="images/plume5c_part_01.png"><img width=200 src="images/plume5c_part_streak_01.png"></a></td>
<td><a href="images/plume5c_part_10.png"><img width=200 src="images/plume5c_part_streak_10.png"></a></td>
<td><a href="images/plume5c_part_30.png"><img width=200 src="images/plume5c_part_streak_30.png"></a></td>
</tr>
</table>

<p><hr>
<a name="3dsmoke"></a>
<p><table>
<tr><th colspan=3> 3D Smoke</th></tr>

<tr>
<td><a href="images/plume5c_part_01.png"><img width=200 src="images/plume5c_smoke_01.png"></a></td>
<td><a href="images/plume5c_part_10.png"><img width=200 src="images/plume5c_smoke_10.png"></a></td>
<td><a href="images/plume5c_part_30.png"><img width=200 src="images/plume5c_smoke_30.png"></a></td>
</tr>
</table>

<p><table>

<tr><td></td><th>all planes</th><th>every 2nd plane</th><th>every 3rd plane</th></tr>
<tr>
<th>0.1 m grid</th>
<td><a href="images/smoke_test_all.png"><img width=300 src="images/smoke_test_all.png"></a></td>
<td><a href="images/smoke_test_every2.png"><img width=300 src="images/smoke_test_every2.png"></a></td>
<td><a href="images/smoke_test_every3.png"><img width=300 src="images/smoke_test_every3.png"></a></td>
</tr>

<tr>
<th> 0.2 m grid</th>
<td><a href="images/smoke_test2_all.png"><img width=300 src="images/smoke_test2_all.png"></a></td>
<td><a href="images/smoke_test2_every2.png"><img width=300 src="images/smoke_test2_every2.png"></a></td>
<td><a href="images/smoke_test2_every3.png"><img width=300 src="images/smoke_test2_every3.png"></a></td>
</tr>
<tr>
<th>actual</th>
<td><a href="images/graysquares.png"><img width=300 src="images/graysquares.png"></a></td>
<td><a href="images/graysquares.png"><img width=300 src="images/graysquares.png"></a></td>
<td><a href="images/graysquares.png"><img width=300 src="images/graysquares.png"></a></td>
</tr>
</table>

<p><hr>
<a name="plot3d"></a>
<p><table>
<tr><th></th><th colspan=4> PLOT3D</th></tr>

<tr>
<th>isosurface</th>
<td><a href="images/plume5c_plot3d_i1.png"><img width=200 src="images/plume5c_plot3d_i1.png"></a></td>
<td><a href="images/plume5c_plot3d_i2.png"><img width=200 src="images/plume5c_plot3d_i2.png"></a></td>
<td><a href="images/plume5c_plot3d_i3.png"><img width=200 src="images/plume5c_plot3d_i3.png"></a></td>
<td><a href="images/plume5c_plot3d_i4.png"><img width=200 src="images/plume5c_plot3d_i4.png"></a></td>
</tr>

<tr>
<th>lines</th>
<td><a href="images/plume5c_plot3d_l1.png"><img width=200 src="images/plume5c_plot3d_l1.png"></a></td>
<td><a href="images/plume5c_plot3d_l2.png"><img width=200 src="images/plume5c_plot3d_l2.png"></a></td>
<td><a href="images/plume5c_plot3d_l3.png"><img width=200 src="images/plume5c_plot3d_l3.png"></a></td>
<td><a href="images/plume5c_plot3d_l4.png"><img width=200 src="images/plume5c_plot3d_l4.png"></a></td>
</tr>

<tr>
<th>vectors</th>
<td><a href="images/plume5c_plot3d_v1.png"><img width=200 src="images/plume5c_plot3d_v1.png"></a></td>
<td><a href="images/plume5c_plot3d_v2.png"><img width=200 src="images/plume5c_plot3d_v2.png"></a></td>
<td><a href="images/plume5c_plot3d_v3.png"><img width=200 src="images/plume5c_plot3d_v3.png"></a></td>
<td><a href="images/plume5c_plot3d_v4.png"><img width=200 src="images/plume5c_plot3d_v4.png"></a></td>
</tr>

<tr><td></td><th>continuous</th><th>stepped</th><th>line</th></tr>
<tr>
<th>shade types</th>
<td><a href="images/plume5c_plot3d_shaded.png"><img width=200 src="images/plume5c_plot3d_shaded.png"></a></td>
<td><a href="images/plume5c_plot3d_step.png"><img width=200 src="images/plume5c_plot3d_step.png"></a></td>
<td><a href="images/plume5c_plot3d_line.png"><img width=200 src="images/plume5c_plot3d_line.png"></a></td>
</tr>
</table>

<p><hr>
<a name="zone"></a>
<p><table>
<tr><th></th><th colspan=2> Zone fire</th></tr>
<tr>
<th>slice</th>
<td><a href="images/cfast_test_c2_200.png"><img width=200 src="images/cfast_test_c2_200.png"></a></d>
<td><a href="images/cfast_test_c2_400.png"><img width=200 src="images/cfast_test_c2_400.png"></a></td>
</tr>

<tr>
<th>3d smoke</th>
<td><a href="images/cfast_test_smoke_200.png"><img width=200 src="images/cfast_test_smoke_200.png"></a></td>
<td><a href="images/cfast_test_smoke_400.png"><img width=200 src="images/cfast_test_smoke_400.png"></a></td>
</tr>
</table>


<p><hr>
<a name="wui"></a>
<p><table>
<tr><th></th><th colspan=4> WUI</th></tr>
<tr>
<th>levelset</th>
<td><a href="images/levelset_30.png"><img width=200 src="images/levelset_30.png"></a></td>
<td><a href="images/levelset_60.png"><img width=200 src="images/levelset_60.png"></a></td>
<td><a href="images/levelset_90.png"><img width=200 src="images/levelset_90.png"></a></td>
<td><a href="images/levelset_120.png"><img width=200 src="images/levelset_120.png"></a></td>
</tr>

<tr>
<th>levelset<br>terrain</th>
<td><a href="images/BT10m_2x2km_LS_0000.png"><img width=200 src="images/BT10m_2x2km_LS_0000.png"></a></td>
<td><a href="images/BT10m_2x2km_LS_0600.png"><img width=200 src="images/BT10m_2x2km_LS_0600.png"></a></td>
<td><a href="images/BT10m_2x2km_LS_1200.png"><img width=200 src="images/BT10m_2x2km_LS_1200.png"></a></td>
<td><a href="images/BT10m_2x2km_LS_1800.png"><img width=200 src="images/BT10m_2x2km_LS_1800.png"></a></td>
</tr>

<tr>
<th>hill</th>
<td><a href="images/hill_structure_015.png"><img width=200 src="images/hill_structure_015.png"></a></td>
<td><a href="images/hill_structure_030.png"><img width=200 src="images/hill_structure_030.png"></a></td>
<td><a href="images/hill_structure_045.png"><img width=200 src="images/hill_structure_045.png"></a></td>
<td><a href="images/hill_structure_060.png"><img width=200 src="images/hill_structure_060.png"></a></td>
</tr>

<tr>
<th>particles </th>
<td><a href="images/pine_tree_part_000.png"><img width=200 src="images/pine_tree_part_000.png"></a></td>
<td><a href="images/pine_tree_part_020.png"><img width=200 src="images/pine_tree_part_020.png"></a></td>
<td><a href="images/pine_tree_part_040.png"><img width=200 src="images/pine_tree_part_040.png"></a></td>
</tr>

<tr>
<th>isosurfaces</th>
<td><a href="images/pine_tree_partiso_000.png"><img width=200 src="images/pine_tree_partiso_000.png"></a></td>
<td><a href="images/pine_tree_partiso_020.png"><img width=200 src="images/pine_tree_partiso_020.png"></a></td>
<td><a href="images/pine_tree_partiso_040.png"><img width=200 src="images/pine_tree_partiso_040.png"></a></td>
</tr>

<tr>
<th>objects</th>
<td><a href="images/tree_test2_part_00.png"><img width=200 src="images/tree_test2_part_00.png"></a></td>
<td><a href="images/tree_test2_part_03.png"><img width=200 src="images/tree_test2_part_03.png"></a></td>
<td><a href="images/tree_test2_part_06.png"><img width=200 src="images/tree_test2_part_06.png"></a></td>
<td><a href="images/tree_test2_part_09.png"><img width=200 src="images/tree_test2_part_09.png"></a></td>
</tr>
<tr>
<th>wind</th>
<td colspan=2><a href="images/wind_test1_002.png"><img width=400 src="images/wind_test1_002.png"></a></td>
</tr>
</table>

<p><hr>
<a name="stereo"></a>
<p><table>
<tr><th></th><th colspan=2>Stereo</th></tr>
<tr>
<th>left/right</th>
<td><img width=375 src="images/thouse5_iso_lr_stereo_L.png"></td>
<td><img width=375 src="images/thouse5_iso_lr_stereo_R.png"></td>
</tr>
<tr>
<th>red/blue</th>
<td colspan=2><img width=450 src="images/thouse5_iso_rb_stereo.png"></td>
</tr>

</table>

<p><hr>
<a name="obstacles"></a>
<p><table>
<tr><th colspan=3> Obstacles</th></tr>
<tr align=center><th>hidden</th><th>outline</th><th>solid</th></tr>
<tr>
<td><a href="images/plume5c_hidden.png"><img width=200 src="images/plume5c_hidden.png"></a></td>
<td><a href="images/plume5c_outline.png"><img width=200 src="images/plume5c_outline.png"></a></td>
<td><a href="images/plume5c_solid.png"><img width=200 src="images/plume5c_solid.png"></a></td>
</tr>
<tr><th colspan=3> Semi-transparent</th></tr>
<tr>
<td><a href="images/transparency_left.png"><img width=200 src="images/transparency_left.png"></a></td>
<td><a href="images/transparency_center.png"><img width=200 src="images/transparency_center.png"></a></td>
<td><a href="images/transparency_right.png"><img width=200 src="images/transparency_right.png"></a></td>
</tr>

<p><hr>
<a name="vents"></a>
<p><table>
<tr><th colspan=3> Rectangular vents</th></tr>
<tr align=center><th>no vents</th><th>no open vents</th><th>all vents</th></tr>
<tr>
<td><a href="images/plume5c_novents.png"><img width=200 src="images/plume5c_novents.png"></a></td>
<td><a href="images/plume5c_noopen.png"><img width=200 src="images/plume5c_noopen.png"></a></td>
<td><a href="images/plume5c_allvents.png"><img width=200 src="images/plume5c_allvents.png"></a></td>
</tr>
</table>

<p><hr>
<p><table>
<tr><th></th><th colspan=4> Circular vents</th></tr>
<tr><td></td><th colspan=2>FDS implementation</th><th colspan=2>user specification</th></tr>
<tr>
<th>solid</th>
<td><a href="images/vcirctest_fds.png"><img width=200 src="images/vcirctest_fds.png"></a></td>
<td><a href="images/vcirctest2_fds.png"><img width=200 src="images/vcirctest2_fds.png"></a></td>
<td><a href="images/vcirctest_circ.png"><img width=200 src="images/vcirctest_circ.png"></a></td>
<td><a href="images/vcirctest2_circ.png"><img width=200 src="images/vcirctest2_circ.png"></a></td>
</tr>
<tr>
<th>outline</th>
<td><a href="images/vcirctest_fds_outline.png"><img width=200 src="images/vcirctest_fds_outline.png"></a></td>
<td><a href="images/vcirctest2_fds_outline.png"><img width=200 src="images/vcirctest2_fds_outline.png"></a></td>
<td><a href="images/vcirctest_circ_outline.png"><img width=200 src="images/vcirctest_circ_outline.png"></a></td>
<td><a href="images/vcirctest2_circ_outline.png"><img width=200 src="images/vcirctest2_circ_outline.png"></a></td>
</tr>
</table>

<p><hr>
<a name="misc"></a>
<p><table>
<tr><th> Texture mapping</th></tr>
<tr><td><a href="images/sillytexture.png"><img width=300 src="images/sillytexture.png"></a></td></tr>
</table>

<p><table>
<tr><th colspan=3> Color conversion</th></tr>
<tr>
<td><a href="images/colorbar_low.png"><img width=300 src="images/colorbar_low.png"></a></td>
<td><a href="images/colorbar_med.png"><img width=300 src="images/colorbar_med.png"></a></td>
<td><a href="images/colorbar_high.png"><img width=300 src="images/colorbar_high.png"></a></td>
</tr>

</table>

<p><table>
<tr><th colspan=3> Clipping</th></tr>
<tr><th>no clipping</th><th>clip blockages not data </th><th>clip blockages and data</th></tr>
<tr>
<td><a href="images/thouse5_smoke_noclip.png"><img width=300 src="images/thouse5_smoke_noclip.png"></a></td>
<td><a href="images/thouse5_smoke_clip_blockages.png"><img width=300 src="images/thouse5_smoke_clip_blockages.png"></a></td>
<td><a href="images/thouse5_smoke_clip_blockages_data.png"><img width=300 src="images/thouse5_smoke_clip_blockages_data.png"></a></td>
</tr>

</table>


<p><hr>
<a name="obj2state"></a>
<p><table>
<tr><th colspan=2> Smokeview objects (two states)</th></tr>
<tr align=center><th>inactive</th><th>active</th></tr>
<tr>
<td><a href="images/sprinkler_inact.png"><img width=200 src="images/sprinkler_inact.png"></a></td>
<td><a href="images/sprinkler_act.png"><img width=200 src="images/sprinkler_act.png"></a></td>
</tr>
<tr>
<td><a href="images/sprinkler_pendent_0.png"><img width=200 src="images/sprinkler_pendent_0.png"></a></td>
<td><a href="images/sprinkler_pendent_1.png"><img width=200 src="images/sprinkler_pendent_1.png"></a></td>
</tr>
<tr>
<td><a href="images/nozzle_0.png"><img width=200 src="images/nozzle_0.png"></a></td>
<td><a href="images/nozzle_1.png"><img width=200 src="images/nozzle_1.png"></a></td>
</tr>
<tr>
<td><a href="images/smokedetector_inact.png"><img width=200 src="images/smokedetector_inact.png"></a></td>
<td><a href="images/smokedetector_act.png"><img width=200 src="images/smokedetector_act.png"></a></td>
</tr>
<tr>
<td><a href="images/heatdetector_inact.png"><img width=200 src="images/heatdetector_inact.png"></a></td>
<td><a href="images/heatdetector_act.png"><img width=200 src="images/heatdetector_act.png"></a></td>
</tr>
</table>

<p><hr>
<a name="obj1state"></a>
<table>
<tr><th colspan=3> Smokeview objects (one state)</th></tr>
<tr>
<td><a href="images/sensor.png"><img width=200 src="images/sensor.png"></a></td>
<td><a href="images/target.png"><img width=200 src="images/target.png"></a></td>
<td><a href="images/fan.png"><img width=200 src="images/fan.png"></a></td>
</tr>
<tr>
<td><a href="images/vent1.png"><img width=200 src="images/vent1.png"></a></td>
<td><a href="images/vent2.png"><img width=200 src="images/vent2.png"></a></td>
<td><a href="images/tube.png"><img width=200 src="images/tube.png"></a></td>
</tr>
</table>

<p><hr>
<a name="objelem"></a>
<p><table>
<tr><th colspan=4>Device elements</th></tr>
<tr>
<td><a href="images/drawarcdisk.png"><img width=200 src="images/drawarcdisk.png"></a></td>
<td><a href="images/drawdisk.png"><img width=200 src="images/drawdisk.png"></a></td>
<td><a href="images/drawcone.png"><img width=200 src="images/drawcone.png"></a></td>
<td><a href="images/drawcube.png"><img width=200 src="images/drawcube.png"></a></td>
</tr>
<tr>
<td><a href="images/drawdisk.png"><img width=200 src="images/drawdisk.png"></a></td>
<td><a href="images/drawhexdisk.png"><img width=200 src="images/drawhexdisk.png"></a></td>
<td><a href="images/drawnotchplate.png"><img width=200 src="images/drawnotchplate.png"></a></td>
<td><a href="images/drawnotchplate2.png"><img width=200 src="images/drawnotchplate2.png"></a></td>
</tr>
<tr>
<td><a href="images/drawpolydisk.png"><img width=200 src="images/drawpolydisk.png"></a></td>
<td><a href="images/drawring.png"><img width=200 src="images/drawring.png"></a></td>
<td><a href="images/drawsphere.png"><img width=200 src="images/drawsphere.png"></a></td>
<td><a href="images/drawtrunccone.png"><img width=200 src="images/drawtrunccone.png"></a></td>
</tr>
</table>
<p><hr>


</BODY>
</HTML>
EOF
