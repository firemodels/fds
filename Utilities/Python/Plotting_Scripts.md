## Introduction

A set of functions has been developed to simplify the generation of plots, which are at the heart of making comparisons between the computation results and the experimental data.  The functions are written in Python (developed in Python 3.7).  If you are new to Python, do not fear, we will go step by step through the process of creating plots.

## Installing Python

The first thing you will need is a functioning version of Python with Matplotlib.  Therefore, unless you otherwise have a preference, I recommend using [pip to manage your Python distribution](https://github.com/firemodels/fds/wiki/Python-Setting-Up-Your-Environment).

Once installed you should have a working version of Python.  Confirm this by opening a [command shell](https://github.com/MaCFP/macfp-db/wiki/Command-Shell) and typing:
```
$ python --version
Python 3.12.2
```
Note that you may need to add or modify the path in your startup script (Linux, macOS).


### Limitations

* The [fdsplotlib.plot_to_fig](https://github.com/MaCFP/macfp-db/wiki/Plotting-Scripts#plot_to_figx_datay_datakwargs) utility, at present, only works for linear plots (log plots will be added shortly).

* The general philosophy in setting up the scripts is that you will be comparing computational results against experimental data.  Thus, each line of a the configuration file is set to plot Exp vs. Cmp.  If you wish to include more than one set of experimental data on a given plot, this is possible, but you need to treat the extra sets of data like computational results. That is, pick one data set to be Exp for all the lines of the configuration file (for that particular plot), then extra datasets can be added and customized as Cmp results.

* Note that the directories _expdir_, _cmpdir_, _pltdir_ may be omitted and the paths included in the filenames _Exp_Filename_ and _Cmp_Filename_. This will be necessary, for example, when we want to make comparisons between various results or experimental data series.

* Please _do not_ use commas, `,`,  or quotation marks, `"` or `'`,  for header labels or values in the configuration file.  Given that we are reading a "comma delimited" spreadsheet file, this simplifies parsing of the file.

### fdsplotlib.dataplot(config_filename,\*\*kwargs)

| Parameter          | Type  | Description |
| ------------------ | ----- | ----------- |
| _config_filename_  | String | (**required**) comma-delimited ASCII text file with plot parameters; each row represents a separate line plot (x,y), see [Configuration File Parameters](); to include more than one line on a single plot, on consecutive rows give the same **Plot_Filename**. |

| Keyword Arguments  | Type  | Description |
| ------------------ | ----- | ----------- |
| _institute_  | String | (optional) [`None`] If present, this string is printed at the top-left of the plot axes |
| _revision_   | String | (optional) [`None`] If present, this string is printed at the top-right of the plot axes |
| _configdir_  | String | (optional) [current directory] Configuration file directory |
| _expdir_     | String | (optional) [current directory] Experimental data directory |
| _cmpdir_     | String | (optional) [current directory] Computational results directory |
| _pltdir_     | String | (optional) [current directory] Plots directory, created if does not exist |
| _close_figs_ | Boolean | (optional) [`False`] Set to `True` will close _all_ open figures after each plot is made.  This can be used to suppress the warning that Python gives for having more than 20 figures open at one time. |
| _verbose_    | Boolean | (optional) [`False`] Set to `True` to print current **Plot_Filename** while processing. |
| _plot_range_  | Range | (optional) [`range(10000)`] A Python [range](https://docs.python.org/3.3/library/stdtypes.html?highlight=range#range) type that can be used to slice the dataframe imported by dataplot for a specific range of rows.  The _start_ and _stop_ for `plot_range` should match the rows in your `*.csv` file, where the header is assumed to be row 1. This is not "Pythonic" but it is much clearer (IMHO) when trying to link dataplot with the configuration file.  Example: `plot_range=range(20,30)`.  Note that `plot_range` is much faster than `plot_list` (next) because it does not have to first read and do string comparisons for the whole configuration file. |
| _plot_list_  | List of Strings | (optional) [`'all'`] A list of **Plot_Filename**'s. Use this to process individual rows of the configuration file. Example: `plot_list=['plot1.pdf','plot2.pdf']` |

Example:
```
fdsplotlib.dataplot(config_filename='FDS_verification_dataplot_config.csv',
                    revision='test',
                    expdir='../../../Experimental_Data/',
                    pltdir='./Plots/')
```

### Configuration File Parameters

Here we give a description of each of the column headings of the configuration file.  These headings set the keyword arguments (\*\*kwargs) for [plot_to_fig]() called by [dataplot]() for both the experimental data and computational results.

Note: some column headers are **required**, while others will take on default values if left out of the configuration file.  Err on the side of being explicit if you want to control the parameters of the plots.

Note: quotes are not needed for strings (including filenames) within the configuration file.  For example, a data label will be given as `Exp`, instead of `"Exp"` or `'Exp'`.

Note: you may **comment** out a line of the configuration file using a hashtag `#`.

| Column heading | Description |
| -------------- | ----------- |
| Exp_Filename   | (**required**) Filename of the ASCII *.csv file for the experimental data (excluding path, if **expdir** is used as a kwarg in dataplot, which is recommended to better organize the config file.) |
| Exp_Header_Row | (optional) [`1`] Integer value of the row used for column headers in the experimental data file.  This may vary (usually 1 or 2) depending on the convention used for placement of units. |
| Exp_Data_Row | (optional) [`Exp_Header_Row+1`] Integer value of the row used for start of experimental data values. |
| Exp_x_Col_Name | (**required**) Column header for experimental independent data (x axis). |
| Exp_y_Col_Name | (**required**) Column header for experimental dependent data (y axis). |
| Exp_y_Error_Col_Name | (optional) [`None`] Column header for experimental data vertical (y axis) uncertainty bands |
| Exp_Data_Markevery | (optional) [`1`] If greater than 1, skip data to make plots more legible. |
| Exp_Error_Absolute | (optional) [`0.`] Absolute data uncertainty range, drawn as transparent band; may be combined with Exp_Error_Relative. `Error_Bound = y_data*(1+/-Exp_Error_Relative)+/-Exp_Error_Absolute`. |
| Exp_Error_Relative | (optional) [`0.`] Relative data uncertainty range, drawn as transparent band; may be combined with Exp_Error_Absolute. `Error_Bound = y_data*(1+/-Exp_Error_Relative)+/-Exp_Error_Absolute`. |
| Exp_Data_Label | (optional) [`Exp`] Label used in legend for experimental data (quotes not needed for strings). |
| Exp_Marker_Style | (optional) [`o`] Experimental data [marker style](https://matplotlib.org/3.3.3/api/markers_api.html#module-matplotlib.markers). Examples: [leave blank], `None`, `o`, `^`, `.`, `+`, `*`, etc.  Note: If you add the Exp_Marker_Style column to your configuration file, then the default will change to whatever you put in that column; a blank will be treated as `None`. |
| Exp_Marker_Edge_Color | (optional) [`black`] Experimental data marker edge color; thickness of edge controlled by Exp_Line_Width below. Examples: `black`, `red`, `blue`, `green` |
| Exp_Marker_Fill_Color | (optional) [`black`] Experimental data marker fill color. Examples: `black`, `red`, `blue`, `green` |
| Exp_Marker_Size | (optional) [`6`] Experimental data marker size. Examples: `0`, `1`, ..., `10` |
| Exp_Line_Style | (optional) [`None`] Experimental data line style. Examples: [leave blank], `None`, `-`, `--`, `-.`, `:` |
| Exp_Line_Color | (optional) [`None`] Experimental data line color. Examples: `black`, `red`, `blue`, `green` |
| Exp_Line_Width | (optional) [`1`] Experimental data line width; also controls marker edge thickness. Examples: `0`, `1`, ..., `10`|
| Cmp_Filename   | (**required**) Filename of the ASCII *.csv file for the computational results |
| Cmp_Header_Row | (optional) [`1`] Integer value of the row used for column headers in the computational results file.  This may vary (usually 1 or 2) depending on the convention used for placement of units. |
| Cmp_Data_Row | (optional) [`Cmp_Header_Row+1`] Integer value of the row used for start of computational results values. |
| Cmp_x_Col_Name | (**required**) Column header for computational independent data (x axis). |
| Cmp_y_Col_Name | (**required**) Column header for computational dependent data (y axis). |
| Cmp_Data_Markevery | (optional) [`1`] If greater than 1, skip data to make plots more legible. |
| Cmp_Data_Label | (optional) [`Cmp`] Label used in legend for computational data (quotes not needed for strings). |
| Cmp_Marker_Style | (optional) [`None`] Computational data [marker style](https://matplotlib.org/3.3.3/api/markers_api.html#module-matplotlib.markers). Examples: [leave blank], `None`, `o`, `^`, `.`, `+`, `*`, etc. |
| Cmp_Marker_Edge_Color | (optional) [`None`] Computational data marker edge color; thickness of edge controlled by Cmp_Line_Width below. Examples: `black`, `red`, `blue`, `green` |
| Cmp_Marker_Fill_Color | (optional) [`None`] Computational data fill marker color. Examples: `black`, `red`, `blue`, `green` |
| Cmp_Marker_Size | (optional) [`0`] Computational data marker size. Examples: `0`, `1`, ..., `10` |
| Cmp_Line_Style | (optional) [`-`] Computational data line style. Examples: [leave blank], `None`, `-`, `--`, `-.`, `:` |
| Cmp_Line_Color | (optional) [`black`] Computational data line color. Examples: `black`, `red`, `blue`, `green` |
| Cmp_Line_Width | (optional) [`1`] Computational data line width; also controls marker edge thickness. Examples: `0`, `1`, ..., `10` |
| Plot_Title | (optional) [` `] Title of the plot, upper-left corner |
| Plot_Subtitle | (optional) [` `] Subtitle of the plot, below title |
| Plot_x_Label | (optional) [`x`] x axis label |
| Plot_y_Label | (optional) [`y`] y axis label |
| Plot_x_Min | (optional) [min(x_data)-5%] Minimum value on x axis |
| Plot_x_Max | (optional) [max(x_data)+5%] Maximum value on x axis |
| Plot_x_Tick | (optional) Tick mark spacing on x axis |
| Plot_x_Nticks | (optional) Number of ticks on x axis (instead of **Plot_x_Tick**) |
| Plot_y_Min | (optional) [min(y_data)-5%] Minimum value on y axis |
| Plot_y_Max | (optional) [max(y_data)+5%] Maximum value on y axis |
| Plot_y_Tick | (optional) Tick mark spacing on y axis |
| Plot_y_Nticks | (optional) Number of ticks on y axis (instead of **Plot_y_Tick**) |
| Plot_x_Axis_Exponent_Min | (optional) [-3] Min exponent before use of scientific notation on x axis tick labels |
| Plot_x_Axis_Exponent_Max | (optional) [3] Max exponent before use of scientific notation on x axis tick labels |
| Plot_y_Axis_Exponent_Min | (optional) [-3] Min exponent before use of scientific notation on y axis tick labels |
| Plot_y_Axis_Exponent_Max | (optional) [3] Max exponent before use of scientific notation on y axis tick labels |
| Plot_Left_Adjust | (optional) Float between 0. and 0.95 used to adjust left side of plot if y axis or tick labels are cutoff; suggest about 0.15, if needed |
| Plot_Bottom_Adjust | (optional) Float between 0. and 0.95 used to adjust bottom side of plot if x axis or tick labels are cutoff; suggest about 0.15, if needed |
| Plot_Right_Adjust | (optional) Float between 0. and 0.95 used to adjust right side of plot if when needed for `outside` legend; suggest about 0.68, if needed |
| Plot_Top_Adjust | (optional) Float between 0. and 0.95 used to adjust top side of plot if when needed |
| Plot_Figure_Width | (optional) [`8`] Usually only needed for **Plot_Legend_Location** `outside`; suggest 12 if needed |
| Plot_Figure_Height | (optional) [`6`] |
| Plot_Flip_Axis | (optional) boolean: `True` or [`False`]; set to `True` to transpose x and y axes data and labels |
| Plot_Show_Legend | (optional) boolean: [`True`] or `False` |
| Plot_Legend_Location | (optional) [`best`] See [matplotlib.pyplot.legend](https://matplotlib.org/3.3.3/api/_as_gen/matplotlib.pyplot.legend.html).  Examples: `best`, `upper right`, `upper left`, `lower right`, `lower left`, `right`, `center left`, `center right`, `lower center`, `upper center`, `center`,`outside`. Note: when using `outside` see **Plot_Right_Adjust** and **Plot_Figure_Width** |
| Plot_Filename | (**required**) Plot filename (excluding path, if **pltdir** is used as a kwarg in dataplot, which is recommended to better organize the config file); equivalent to _fname_ in [matplotlib.pyplot.savefig](https://matplotlib.org/3.3.3/api/_as_gen/matplotlib.pyplot.savefig.html); include filename extension. Examples: `myplot.pdf`, `myplot.png` |

### plot_to_fig(x_data,y_data,**kwargs)

This function is a lower level function that uses [Matplotlib Pyplot](https://matplotlib.org/api/_as_gen/matplotlib.pyplot.html) to generate a simple *x,y* plot and then returns the figure handle.  You can then use [plt.show()](https://matplotlib.org/api/_as_gen/matplotlib.pyplot.show.html) to view the figure and/or [plt.save_fig()](https://matplotlib.org/api/_as_gen/matplotlib.pyplot.savefig.html) to save the figure (this is done automatically by [dataplot](https://github.com/MaCFP/macfp-db/wiki/Plotting-Scripts#macfpdataplotconfig_filenamekwargs)).

| Parameter          | Type    | Description |
| ------------------ | ------- | ----------- |
| _x_data_           | Numeric | (required) x data read from csv file using [pandas.read_csv()](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_csv.html) |
| _y_data_           | Numeric | (required) y data (same size as x data) read from csv file |

| Keyword Arguments  | Type  | Description |
| ------------------ | ----- | ----------- |
| _data_markevery_   | Integer | (optional) [`1`] Skip this number of points when plotting.  |
| _data_label_       | String  | (optional) [`None`] Label for data in legend. |
| _y_error_absolute_ | Numeric | (optional) [`0.`] Absolute error |
| _y_error_relative_ | Numeric | (optional) [`0.`] Relative error, may be used in place of absolute error |
| _y_error_vector_ | Numeric | (optional) [`0.`] y data error, must be same size as _y_data_ |
| _x_label_       | String  | (optional) [`None`] Label for x axis  |
| _y_label_       | String  | (optional) [`None`] Label for y axis |
| _marker_style_       | String  | (optional) [`None`] Examples: `o`, `*`, `+`, `>`, `<`, `^`, `v`, `.`, `s` (square), [etc.](https://matplotlib.org/3.3.3/api/markers_api.html#module-matplotlib.markers)  |
| _marker_edge_color_       | String  | (optional) [`black`] Examples: `black`, `red`, `blue`, `green`  |
| _marker_fill_color_       | String  | (optional) [`None`] Examples: `black`, `red`, `blue`, `green`  |
| _marker_size_       | Integer  | (optional) [`6`] `1`, `2`, etc. (`0` not allowed, just leave blank) |
| _line_style_       | String  | (optional) [`None`] Examples: `-`, `--`, `-.`, `:`  |
| _line_color_       | String  | (optional) [`None`] Examples: `black`, `red`, `blue`, `green`  |
| _line_width_       | Integer | (optional) [`1`] `1`, `2`, etc. (`0` not allowed, just leave blank) |
| _x_min_       | Numeric  | (optional) [`None`] Minimum value on x axis |
| _x_max_       | Numeric  | (optional) [`None`] Maximum value on x axis |
| _n_xticks_       | Integer  | (optional) [`None`] Number of ticks on x axis |
| _y_min_       | Numeric  | (optional) [`None`] Minimum value on y axis |
| _y_max_       | Numeric  | (optional) [`None`] Maximum value on y axis |
| _n_yticks_       | Integer  | (optional) [`None`] Number of ticks on y axis |
| _show_legend_       | String  | (optional) [`None`] boolean: True or [False] |
| _legend_location_       | String  | (optional) [`None`] Examples: `best`, `upper right`, `upper left`,`outside` (only outside right available; see _figure_right_adjust_; apply to all calls to _plot_to_fig_ with the same figure handle) |
| _legend_framealpha_       | Numeric  | (optional) [`0.8`] Alpha (opposite of transparency) of legend background. Change to 1.0 for no transpency. |
| _ticklabel_fontsize_       | Integer  | (optional) [`16`], `24`, etc. |
| _axeslabel_fontsize_       | Integer  | (optional) [`18`], `24`, etc. |
| _title_fontsize_       | Integer  | (optional) [`18`], `24`, etc. |
| _subtitle_fontsize_       | Integer  | (optional) [`16`], `24`, etc. |
| _legend_fontsize_       | Integer  | (optional) [`16`], `24`, etc. |
| _plot_title_       | String  | (optional) [`None`] Plot title upper left |
| _plot_subtitle_       | String  | (optional) [`None`] Plot subtitle below title |
| _institute_label_       | String  | (optional) [`None`] Stamp upper left above axes |
| _figure_size_       | Numeric Tuple  | (optional) [`(8,6)`] (Width, Height) in inches |
| _figure_handle_       | Integer  | (optional) plot to a specific figure |
| _figure_left_adjust_  | Numeric  | (optional) 0<x<0.95; use to adjust left if y axis or tick labels are cut off |
| _figure_bottom_adjust_ | Numeric  | (optional) 0<x<0.95; use to adjust bottom if x axis or tick labels are cut off |
| _figure_right_adjust_  | Numeric  | (optional) 0<x<0.95; use to adjust right with using `legend_location='outside'` |
| _figure_top_adjust_   | Numeric  | (optional) 0<x<0.95; use to adjust top |
| _figure_x_axis_exponent_min_   | Integer  | (optional) [-3]; set min threshold for x axis tick labels to use scientific notation |
| _figure_x_axis_exponent_max_   | Integer  | (optional) [3]; set max threshold for x axis tick labels to use scientific notation |
| _figure_y_axis_exponent_min_   | Integer  | (optional) [-3]; set min threshold for y axis tick labels to use scientific notation |
| _figure_y_axis_exponent_max_   | Integer  | (optional) [3]; set max threshold for y axis tick labels to use scientific notation |




