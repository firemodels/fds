#!/usr/bin/env python
"""
fdsplotlib.py
by Randy McDermott
January 2025

Fire Dynamics Simulator (FDS) Plot Library

Collection of functions for plotting and analysis
"""

import os
import sys
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
# rc('text', usetex=True) # Enable TeX rendering
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd


def expand_ranges(items, df, header_rows=1):
    """
    Expand a list of row specifiers into pandas row indices.
    Supports:
      - int: row number (1-based, with headers counted)
      - "start:stop": inclusive range of row numbers
      - "start:": open-ended range (to the end)
      - "all": keep everything
      - plain string: match Dataname column (case-insensitive)
    """
    nrows = len(df)
    result = []

    for item in items:
        if isinstance(item, int):
            # single row number
            result.append(item - (header_rows + 1))

        elif isinstance(item, str):
            s = item.strip()
            low = s.lower()

            if low == "all":
                return list(range(nrows))  # everything

            if ":" in s:  # range form
                start, _, end = s.partition(":")
                start = int(start)

                if end == "":
                    # open-ended: go to the last row of df
                    end = nrows + header_rows  # Excel-style row count
                else:
                    end = int(end)

                # Convert to iloc positions (0-based, exclusive of end)
                start_pos = start - (header_rows + 1)
                end_pos = end - header_rows

                # Clamp so we don't go past the last row
                end_pos = min(end_pos, nrows)

                rng = range(start_pos, end_pos)
                result.extend(rng)

            else:
                # assume it's a Dataname match
                matches = df.index[df['Dataname'].str.lower() == low].tolist()
                if not matches:
                    raise ValueError(f"No match for Dataname '{item}'")
                result.extend(matches)

        else:
            raise TypeError(f"Unsupported plot_range element: {item}")

    return sorted(set(result))


def dataplot(config_filename,**kwargs):

    import os
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib import rc

    import logging
    # Suppress just the 'findfont' warnings from matplotlib's font manager
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

    # defaults

    configdir = kwargs.get('configdir','')
    revision = kwargs.get('revision','')
    expdir = kwargs.get('expdir','')
    cmpdir = kwargs.get('cmpdir','')
    pltdir = kwargs.get('pltdir','')
    close_figs = kwargs.get('close_figs', False)
    verbose = kwargs.get('verbose', False)

    plot_list = kwargs.get('plot_list', ['all'])
    plot_range_in = kwargs.get('plot_range', None)
    header_rows = kwargs.get('header_rows', 1)

    # --- Initialize scaffolding for scatplot compatibility ---
    drange = []
    Save_Quantity = []
    Save_Group_Style = []
    Save_Fill_Color = []
    Save_Group_Key_Label = []
    Save_Measured_Metric = []
    Save_Predicted_Metric = []
    Save_Dataname = []
    Save_Plot_Filename = []
    Save_Dep_Title = []
    Save_Error_Tolerance = []
    Save_Metric_Type = []
    Save_Measured_Quantity = []
    Save_Predicted_Quantity = []

    # Rebuild the default NA strings *excluding* 'N/A'
    default_na = {
        '', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN',
        '1.#IND', '1.#QNAN', '<NA>', 'NA', 'NULL', 'NaN',
        'n/a', 'nan', 'null'
    }
    safe_na_values = default_na  # 'N/A' intentionally not included

    # read the config file
    df = pd.read_csv(configdir+config_filename, sep=',', engine='python', quotechar='"', na_values=safe_na_values, keep_default_na=False)
    C = df.where(pd.notnull(df), None)

    if plot_range_in is not None:
        if isinstance(plot_range_in, (list, tuple)):
            adjusted = expand_ranges(plot_range_in, C, header_rows)
            if len(adjusted) < len(C):  # filtering happened
                C = C.iloc[adjusted]
        elif isinstance(plot_range_in, range):
            start, end = plot_range_in.start, plot_range_in.stop
            adjusted = range(start - (header_rows + 1),
                             end - header_rows)  # Python range already exclusive
            C = C.iloc[adjusted]
        else:
            raise TypeError("plot_range must be a list, tuple, or range")
    elif plot_list and 'all' not in [p.lower() for p in plot_list]:
        # Only filter by plot_list if no plot_range was passed
        C = C[C['Dataname'].str.lower().isin([p.lower() for p in plot_list])]

    Plot_Filename_Last = None
    d1_Key_Last = None
    f_Last = plt.figure()

    # loop over the rows of the config file
    for pos, (irow, row) in enumerate(C.iterrows()):

        pp = define_plot_parameters(C, pos)  # use position, not label

        # debug dump
        # print(vars(pp))

        # print(pp.__dict__) # helpful for debug

        # ----------------------------------------------------------------------
        # Handle MATLAB dataplot switch_id behavior (d, f, o, g, s)
        # ----------------------------------------------------------------------
        switch_id = str(pp.switch_id).strip().lower()

        # Skip 's' outright
        if switch_id == 's':
            continue

        # If ANY 'o' lines exist in the filtered config C, process only those
        otest_active = any( str(C.iloc[j]['switch_id']).strip().lower() == 'o' for j in range(len(C)) )

        if otest_active and switch_id != 'o':
            continue

        # 'g' lines: generate plots but EXCLUDE from scatplot stats & drange
        gtest = (switch_id == 'g')

        # 'd' default, 'f' follow-on (your filename reuse already mimics “hold on”)
        dtest = (switch_id == 'd')
        ftest = (switch_id == 'f')

        # If it’s none of the recognized ones, skip safely
        if not (dtest or ftest or gtest or switch_id == 'o'):
            if verbose:
                print(f"[dataplot] Skipping unrecognized switch_id '{pp.switch_id}' on line {irow+2}")
            continue

        # Track drange like MATLAB (1-based CSV lines starting at row 2)
        if not gtest:
            drange.append(irow + 2)

        # Append metadata only for rows that should appear in scatplot
        if not gtest:
            Save_Dataname.append(pp.Dataname)
            Save_Plot_Filename.append(pp.Plot_Filename)
            Save_Dep_Title.append(pp.Dep_Title)
            Save_Error_Tolerance.append(pp.Error_Tolerance)
            Save_Metric_Type.append(pp.Metric)
            Save_Group_Key_Label.append(pp.Group_Key_Label)
            Save_Quantity.append(pp.Quantity)
            Save_Group_Style.append(pp.Group_Style)
            Save_Fill_Color.append(pp.Fill_Color)

            # Placeholders (we will overwrite for this row below)
            Save_Measured_Metric.append(np.nan)
            Save_Predicted_Metric.append(np.nan)
            Save_Measured_Quantity.append(None)
            Save_Predicted_Quantity.append(None)



        if pp.Plot_Filename!=Plot_Filename_Last:

            if verbose:
                print(f'Generating plot {irow+2} ' + pltdir + pp.Plot_Filename + '...')

            if close_figs:
                plt.close('all')

            # read data from exp file
            # set header to the row where column names are stored (Python is 0 based)
            E = pd.read_csv(expdir+pp.d1_Filename, header=int(pp.d1_Col_Name_Row-1), sep=',', engine='python', quotechar='"')
            E.columns = E.columns.str.strip()  # <-- Strip whitespace from headers

            start_idx = int(pp.d1_Data_Row - pp.d1_Col_Name_Row - 1)
            x, col_names = get_data(E, pp.d1_Ind_Col_Name, start_idx)
            y, col_names = get_data(E, pp.d1_Dep_Col_Name, start_idx)

            if pp.d1_Style:
                raw_styles = [c.strip() for c in pp.d1_Style.split('|')]
            else:
                raw_styles = []
            styles = (raw_styles + [None] * len(col_names))[:len(col_names)]

            if pp.d1_Key:
                raw_keys = [c.strip() for c in pp.d1_Key.split('|')]
            else:
                raw_keys = []
            # Pad or truncate to match col_names
            key_labels = (raw_keys + [None] * len(col_names))[:len(col_names)]

            for i, label in enumerate(col_names):
                if i==0:
                    # plot the exp data
                    f = plot_to_fig(x_data=x, y_data=y[:, i],
                        data_label=key_labels[i],
                        x_label=pp.Ind_Title,
                        y_label=pp.Dep_Title,
                        marker_style=styles[i],
                        x_min=pp.Min_Ind,x_max=pp.Max_Ind,
                        y_min=pp.Min_Dep,y_max=pp.Max_Dep,
                        legend_location=matlab_legend_to_matplotlib(pp.Key_Position)
                        )
                elif i>0:
                    f = plot_to_fig(x_data=x, y_data=y[:, i],
                        figure_handle=f,
                        data_label=key_labels[i],
                        x_label=pp.Ind_Title,
                        y_label=pp.Dep_Title,
                        marker_style=styles[i],
                        x_min=pp.Min_Ind,x_max=pp.Max_Ind,
                        y_min=pp.Min_Dep,y_max=pp.Max_Dep,
                        legend_location=matlab_legend_to_matplotlib(pp.Key_Position)
                        )

            # plt.figure(f.number) # make figure current
            # plt.show()
        else:
            f = f_Last

            if pp.d1_Key!=d1_Key_Last:

                # read data from exp file
                # set header to the row where column names are stored (Python is 0 based)
                E = pd.read_csv(expdir+pp.d1_Filename, header=int(pp.d1_Col_Name_Row-1), sep=',', engine='python', quotechar='"')
                E.columns = E.columns.str.strip()  # <-- Strip whitespace from headers
                start_idx = int(pp.d1_Data_Row - pp.d1_Col_Name_Row - 1)
                x, col_names = get_data(E, pp.d1_Ind_Col_Name, start_idx)
                y, col_names = get_data(E, pp.d1_Dep_Col_Name, start_idx)

                # plot the exp data
                f = plot_to_fig(x_data=x, y_data=y,
                    figure_handle=f,
                    data_label=pp.d1_Key,
                    x_label=pp.Ind_Title,
                    y_label=pp.Dep_Title,
                    marker_style=pp.d1_Style,
                    x_min=pp.Min_Ind,x_max=pp.Max_Ind,
                    y_min=pp.Min_Dep,y_max=pp.Max_Dep,
                    legend_location=matlab_legend_to_matplotlib(pp.Key_Position)
                    )

        # --- Save measured (experimental) metric using the numeric arrays from get_data ---
        try:
            if not gtest:
                Save_Measured_Metric[-1] = compute_metric(y, pp.Metric, x, pp.d1_Initial_Value)
                Save_Measured_Quantity[-1] = col_names
        except Exception as e:
            print(f"[dataplot] Warning: measured metric failed for {pp.Dataname}: {e}")


        # get the model results
        M = pd.read_csv(cmpdir+pp.d2_Filename, header=int(pp.d2_Col_Name_Row-1), sep=',', engine='python', quotechar='"')
        M.columns = M.columns.str.strip()  # <-- Strip whitespace from headers
        start_idx = int(pp.d2_Data_Row - pp.d2_Col_Name_Row - 1)

        x, col_names = get_data(M, pp.d2_Ind_Col_Name, start_idx)
        y, col_names = get_data(M, pp.d2_Dep_Col_Name, start_idx)

        version_string = revision
        if (pp.VerStr_Filename):
            file1 = open(cmpdir+pp.VerStr_Filename,"r")
            Lines = file1.readlines()
            version_string = Lines[0].strip()
            file1.close()

        if pp.d2_Style:
            raw_styles = [c.strip() for c in pp.d2_Style.split('|')]
        else:
            raw_styles = []
        styles = (raw_styles + [None] * len(col_names))[:len(col_names)]

        if pp.d2_Key:
            raw_keys = [c.strip() for c in pp.d2_Key.split('|')]
        else:
            raw_keys = []
        key_labels = (raw_keys + [None] * len(col_names))[:len(col_names)]

        for i, label in enumerate(col_names):
            f = plot_to_fig(x_data=x, y_data=y[:, i],
                revision_label=version_string,
                figure_handle=f,
                x_label=pp.Ind_Title,
                y_label=pp.Dep_Title,
                data_label=key_labels[i],
                line_style=styles[i],
                x_min=pp.Min_Ind,x_max=pp.Max_Ind,
                y_min=pp.Min_Dep,y_max=pp.Max_Dep,
                legend_location=matlab_legend_to_matplotlib(pp.Key_Position),
                plot_title=pp.Plot_Title
                )


        # --- Save predicted (model) metric using the numeric arrays from get_data ---
        try:
            if not gtest:
                Save_Predicted_Metric[-1] = compute_metric(y, pp.Metric, x, pp.d2_Initial_Value)
                Save_Predicted_Quantity[-1] = col_names
        except Exception as e:
            print(f"[dataplot] Warning: predicted metric failed for {pp.Dataname}: {e}")



        plt.figure(f.number) # make figure current
        # plt.show()

        # create plot directory if it does not exist
        isDir = os.path.isdir(pltdir)
        if not isDir:
            os.mkdir(pltdir)

        plt.savefig(pltdir + pp.Plot_Filename + '.pdf', backend='pdf')

        Plot_Filename_Last = pp.Plot_Filename
        d1_Key_Last = pp.d1_Key
        f_Last = f

        # except:
        #     print("Error in row {whichrow}, skipping case...".format(whichrow=irow+1))
        #     continue

    # --- MATLAB-compatible output scaffolding for scatplot interface ---
    try:
        saved_data = [
            Save_Quantity,
            Save_Group_Style,
            Save_Fill_Color,
            Save_Group_Key_Label,
            Save_Measured_Metric,
            Save_Predicted_Metric,
            Save_Dataname,
            Save_Plot_Filename,
            Save_Dep_Title,
            Save_Error_Tolerance,
            Save_Metric_Type,
            Save_Measured_Quantity,
            Save_Predicted_Quantity,
        ]
    except Exception as e:
        print(f"[dataplot] Error assembling saved_data: {e}")
        saved_data = []

    print("[dataplot] returning saved_data and drange")
    return saved_data, drange


def compute_metric(data, metric_type, x=None, initial_value=0.0):
    import numpy as np
    try:
        if metric_type == 'max':
            return np.nanmax(data) - initial_value
        elif metric_type == 'min':
            return initial_value - np.nanmin(data)
        elif metric_type == 'maxabs':
            return np.nanmax(np.abs(data - initial_value))
        elif metric_type == 'mean':
            return abs(np.nanmean(data) - initial_value)
        elif metric_type == 'area' and x is not None:
            return np.trapz(data, x) - initial_value
        elif metric_type == 'all':
            return data - initial_value
        else:
            return np.nanmean(data) - initial_value  # default
    except Exception:
        return np.nan


def get_data(E, spec, start_idx):
    """
    Extract data columns from DataFrame E according to spec string.

    spec: "colA" | "colA|colB" | "colA+colB"
    start_idx: integer row index where numeric data starts

    Returns
    -------
    tuple (y, col_names)
        y : 2D numpy array (nrows x ncols)
        col_names : list of str, resolved column names (with '+' grouped)
    """
    names = [s.strip() for s in spec.split('|')]
    out = []
    col_names = []
    for name in names:
        if '+' in name:
            cols = [n.strip() for n in name.split('+')]
            series = E[cols].iloc[start_idx:].astype(float).sum(axis=1).values
            out.append(series)
            col_names.append('+'.join(cols))  # keep group name
        else:
            series = E[[name]].iloc[start_idx:].astype(float).values.ravel()
            out.append(series)
            col_names.append(name)
    y = np.column_stack(out) if len(out) > 1 else np.array(out[0]).reshape(-1, 1)
    return y, col_names


def plot_to_fig(x_data,y_data,**kwargs):
    """
    Create a simple x,y plot and return the fig handle
    """
    # # useful debug statements
    # print(x_data)
    # print(y_data)
    # for key, value in kwargs.items():
    #     print ("%s == %s" %(key, value))

    plot_style = get_plot_style("fds")
    plt.rcParams["font.family"] = plot_style["Font_Name"]
    plt.rcParams["font.size"] = plot_style["Label_Font_Size"]
    # print(plot_style)

    plt.rcParams['text.usetex'] = True # supports latex math (set per plot below)
    plt.rcParams["pdf.use14corefonts"] = True # forces matplotlib to write native pdf fonts rather than embed
    import logging
    # Suppress just the 'findfont' warnings from matplotlib's font manager
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

    ##### default parameters ######
    default_figure_size = (plot_style["Paper_Width"],plot_style["Paper_Height"])
    default_plot_size = (plot_style["Plot_Width"],plot_style["Plot_Height"])
    default_plot_origin = (plot_style["Plot_X"],plot_style["Plot_Y"])
    default_ticklabel_fontsize = plot_style["Label_Font_Size"]
    default_axeslabel_fontsize = plot_style["Label_Font_Size"]
    default_legend_fontsize = plot_style["Key_Font_Size"]
    default_title_fontsize = plot_style["Title_Font_Size"]
    default_version_fontsize = 10
    default_legend_location = 'best'
    default_legend_framealpha = 1
    default_markevery = 1
    markerfacecolor = None
    markeredgecolor = 'black'
    markeredgewidth = 1
    marker = None
    markersize = 5
    linestyle = '-'
    linewidth = 1
    color = 'black'
    ###############################

    figure_size=kwargs.get('figure_size',default_figure_size)
    plot_size=kwargs.get('plot_size',default_plot_size)
    plot_origin=kwargs.get('plot_origin',default_plot_origin)
    version_fontsize=kwargs.get('version_fontsize',default_version_fontsize)

    # if figure handle is passed, append to current figure, else generate a new figure
    if kwargs.get('figure_handle'):
        fig = kwargs.get('figure_handle')
        if fig.axes:
            ax = fig.axes[0]
        else:
            # Ensure we have at least one axes
            ax = fig.add_subplot(111)
        plt.figure(fig.number)
        using_existing_figure = True
    else:
        fig = plt.figure(figsize=figure_size)
        using_existing_figure = False
        # Convert to fractions of the figure size:
        ax_w = plot_size[0] / figure_size[0]
        ax_h = plot_size[1] / figure_size[1]
        left   = plot_origin[0] / figure_size[0]
        bottom = plot_origin[1] / figure_size[1]
        ax = fig.add_axes([left, bottom, ax_w, ax_h])

    # select plot type
    plot_type=kwargs.get('plot_type','linear')

    # convert matlab styles to matplotlib
    style = kwargs.get('marker_style','ko')
    color,marker,linestyle = parse_matlab_style(style)

    if kwargs.get('line_style'):
        style = kwargs.get('line_style')
        color,marker,linestyle = parse_matlab_style(style)

    marker_fill_color = kwargs.get('marker_fill_color',None)
    markerfacecolor = marker_fill_color

    error_fill_color = kwargs.get('error_fill_color',None)

    # other plot parameters
    markevery = kwargs.get('data_markevery',default_markevery)
    legend_location = kwargs.get('legend_location',default_legend_location)
    legend_framealpha = kwargs.get('legend_framealpha',default_legend_framealpha)

    data_label = kwargs.get('data_label',None)

    # trap any data_labels set to blank (old matlab convention)
    if isinstance(data_label, str) and data_label.lower() == 'blank':
        data_label = None

    # generate the main x,y plot
    if plot_type=='linear':
        ax.plot(x_data,y_data,
            markevery=markevery,
            label=data_label,
            markerfacecolor=markerfacecolor,
            markeredgecolor=color,
            markeredgewidth=markeredgewidth,
            marker=marker,
            markersize=markersize,
            linestyle=linestyle,
            linewidth=linewidth,
            color=color)

    if plot_type=='loglog':
        ax.loglog(x_data,y_data,
            markevery=markevery,
            label=data_label,
            markerfacecolor=markerfacecolor,
            markeredgecolor=color,
            markeredgewidth=markeredgewidth,
            marker=marker,
            markersize=markersize,
            linestyle=linestyle,
            linewidth=linewidth,
            color=color)

    if plot_type=='semilogx':
        ax.semilogx(x_data,y_data,
            markevery=markevery,
            label=data_label,
            markerfacecolor=markerfacecolor,
            markeredgecolor=color,
            markeredgewidth=markeredgewidth,
            marker=marker,
            markersize=markersize,
            linestyle=linestyle,
            linewidth=linewidth,
            color=color)

    if plot_type=='semilogy':
        ax.semilogy(x_data,y_data,
            markevery=markevery,
            label=data_label,
            markerfacecolor=markerfacecolor,
            markeredgecolor=color,
            markeredgewidth=markeredgewidth,
            marker=marker,
            markersize=markersize,
            linestyle=linestyle,
            linewidth=linewidth,
            color=color)

    # if y fill range is passed, add it to the plot
    if kwargs.get('y_error_fill_absolute') and not kwargs.get('y_error_fill_relative'):
        if kwargs.get('y_error_fill_absolute')>0.:
            ax.fill_between(x_data,y_data-kwargs.get('y_fill_absolute'),y_data+kwargs.get('y_fill_absolute'),
                alpha=0.1,color=error_fill_color)

    if kwargs.get('y_error_fill_relative') and not kwargs.get('y_error_fill_absolute'):
        if kwargs.get('y_error_fill_relative')>0.:
            ax.fill_between(x_data,y_data*(1.-kwargs.get('y_error_fill_relative')),y_data*(1.+kwargs.get('y_error_fill_relative')),
                alpha=0.1,color=error_fill_color)

    if kwargs.get('y_error_fill_relative') and kwargs.get('y_error_fill_absolute'):
        if kwargs.get('y_error_fill_relative')>0.:
            ax.fill_between(x_data,y_data*(1.-kwargs.get('y_error_fill_relative'))-kwargs.get('y_error_fill_absolute'),
                                   y_data*(1.+kwargs.get('y_error_fill_relative'))+kwargs.get('y_error_fill_absolute'),
                alpha=0.1,color=error_fill_color)

    if kwargs.get('y_error_fill'):
        y_error_fill = kwargs.get('y_error_fill')
        if len(y_data)==len(y_error_fill):
            ax.fill_between(x_data,y_data-y_error_fill,y_data+y_error_fill,
                alpha=0.1,color=error_fill_color)
        else:
            raise ValueError(f"y_fill must the same length as y_data")

    xerr = kwargs.get('x_error', None)
    yerr = kwargs.get('y_error', None)
    if xerr is not None or yerr is not None:
        ax.errorbar(
            x_data, y_data,
            xerr=xerr,                               # can be scalar, array, or [lower, upper]
            yerr=yerr,                               # same flexibility
            fmt=style,                               # marker style for data points
            markeredgewidth=markeredgewidth,         # marker edge width
            markerfacecolor=markerfacecolor,         # make marker hollow
            markeredgecolor=color,                   # outline color
            linestyle=linestyle,
            linewidth=linewidth,
            capsize=kwargs.get('error_capsize', 5),  # size of caps at ends
            capthick=linewidth,
        )


    ticklabel_fontsize=kwargs.get('ticklabel_fontsize',default_ticklabel_fontsize)
    plt.setp( ax.xaxis.get_majorticklabels(), rotation=0, fontsize=ticklabel_fontsize )
    plt.setp( ax.yaxis.get_majorticklabels(), rotation=0, fontsize=ticklabel_fontsize )

    axeslabel_fontsize=kwargs.get('axeslabel_fontsize',default_axeslabel_fontsize)
    if not using_existing_figure:
        plt.xlabel(kwargs.get('x_label'), fontsize=axeslabel_fontsize)
        plt.ylabel(kwargs.get('y_label'), fontsize=axeslabel_fontsize)

    legend_fontsize=kwargs.get('legend_fontsize',default_legend_fontsize)

    if data_label:
        if kwargs.get('legend_location')=='outside':
            plt.legend(fontsize=legend_fontsize,bbox_to_anchor=(1,1),loc='upper left',framealpha=legend_framealpha)
        else:
            # if kwargs.get('show_legend'):
            plt.legend(fontsize=legend_fontsize,loc=legend_location,framealpha=legend_framealpha)

    # plot title
    if kwargs.get('plot_title'):
        if kwargs.get('title_fontsize'):
            title_fontsize=kwargs.get('title_fontsize')
        else:
            title_fontsize=default_title_fontsize

        plt.text(0.05, 0.95, kwargs.get('plot_title'),
        transform=plt.gca().transAxes,
        fontsize=title_fontsize,
        verticalalignment='top',
        horizontalalignment='left')

    # set axes and tick properties
    xmin=kwargs.get('x_min')
    xmax=kwargs.get('x_max')
    ymin=kwargs.get('y_min')
    ymax=kwargs.get('y_max')

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    if plot_type in ('loglog', 'semilogy'):
        apply_global_exponent(ax, axis='y', fontsize=axeslabel_fontsize)
    if plot_type in ('loglog', 'semilogx'):
        apply_global_exponent(ax, axis='x', fontsize=axeslabel_fontsize)

    if kwargs.get('revision_label'):
        add_version_string(ax=ax, version_str=kwargs.get('revision_label'), plot_type=plot_type, font_size=version_fontsize)

    # fig.tight_layout() # this should not be needed if figure_size and plot_size are both specified

    return fig


def apply_global_exponent(ax, axis='y', fontsize=10, minor_subs=None):
    import numpy as np
    from matplotlib.ticker import FuncFormatter, LogLocator

    if axis == 'y':
        ticks = ax.get_yticks()
        axis_obj = ax.yaxis
        lims = ax.get_ylim()
    else:
        ticks = ax.get_xticks()
        axis_obj = ax.xaxis
        lims = ax.get_xlim()

    # Keep only positive finite ticks (for log axes)
    ticks = np.array([t for t in ticks if t > 0 and np.isfinite(t)])
    if ticks.size == 0:
        return

    # Choose representative exponent
    exp = int(np.floor(np.log10(np.median(ticks))))
    scale = 10.0**exp

    # Major formatter: fixed-point decimals
    def fmt(val, pos):
        v = val / scale
        return "{:g}".format(v)

    axis_obj.set_major_formatter(FuncFormatter(fmt))
    axis_obj.get_offset_text().set_visible(False)

    # Decide what to do with minor tick labels
    span_decades = np.log10(lims[1]) - np.log10(max(lims[0], 1e-300))

    # Default subs = 2, 4, 6, 8
    if minor_subs is None:
        minor_subs = [2, 4, 6, 8]

    if span_decades <= 1.1:  # only ~1 decade
        axis_obj.set_minor_locator(LogLocator(base=10.0, subs=minor_subs, numticks=10))
        axis_obj.set_minor_formatter(FuncFormatter(lambda val, pos: "{:g}".format(val/scale)))
    else:
        ax.tick_params(axis=axis, which='minor', labelleft=False, labelbottom=False)

    # Force tick labels NOT to go through TeX
    if axis == 'y':
        for label in ax.get_yticklabels() + ax.get_yticklabels(minor=True):
            label.set_usetex(False)
    else:
        for label in ax.get_xticklabels() + ax.get_xticklabels(minor=True):
            label.set_usetex(False)

    # Place the ×10^exp text at the axis end
    if exp != 0:
        if axis == 'y':
            ax.text(0, 1.01, rf"$\times 10^{{{exp}}}$",
                    transform=ax.transAxes,
                    ha='left', va='bottom', fontsize=fontsize)
        else:
            ax.text(1.0, -0.1, rf"$\times 10^{{{exp}}}$",
                    transform=ax.transAxes,
                    ha='right', va='top', fontsize=fontsize)



def parse_matlab_style(style):
    color = ''
    marker = ''
    linestyle = ''

    # Check for the color (first character)
    if style[0] == 'k':
        color = 'black'
    elif style[0] == 'r':
        color = 'red'
    elif style[0] == 'g':
        color = 'green'
    elif style[0] == 'b':
        color = 'blue'
    elif style[0] == 'y':
        color = 'yellow'
    elif style[0] == 'm':
        color = 'magenta'
    elif style[0] == 'c':
        color = 'cyan'
    elif style[0] == 'w':
        color = 'white'
    else:
        raise ValueError(f"Unknown color code: {style[0]}")

    # Check for the marker style (rest of the string)
    for char in style[1:]:
        if char == 'o':
            marker = 'o'  # Circle
        elif char == 's':
            marker = 's'  # Square
        elif char == 'd':
            marker = 'd'  # Diamond
        elif char == '^':
            marker = '^'  # Triangle up
        elif char == 'v':
            marker = 'v'  # Triangle down
        elif char == '>':
            marker = '>'  # Triangle right
        elif char == '<':
            marker = '<'  # Triangle left
        elif char == '*':
            marker = '*'  # Star
        elif char == '+':
            marker = '+'  # Plus
        elif char == 'x':
            marker = 'x'  # X
        # Ignore unknown characters (to handle linestyle separately)

    # Check for the line style
    if '--' in style:
        linestyle = '--'  # Dashed
    elif '-.' in style:
        linestyle = '-.'  # Dash-dot
    elif '-' in style:
        linestyle = '-'   # Solid
    elif ':' in style:
        linestyle = ':'   # Dotted
    else:
        linestyle = ''    # No line style

    return color, marker, linestyle


def get_version_string(filename):
    file1 = open(filename,"r")
    Lines = file1.readlines()
    version_str = Lines[0].strip()
    file1.close()
    return version_str

def add_version_string(ax, version_str, plot_type='linear', scale_x=1.00, scale_y=1.02,
                       font_name='Times', font_size=10):
    """
    Adds a version string to a matplotlib plot.

    Parameters:
    ax (matplotlib.axes.Axes): The axes to add the version string to.
    filename (str): Path to the version string file.
    plot_type (str): Type of plot ('loglog', 'semilogx', 'semilogy', or 'linear').
    scale_x (float): Scaling factor for X-position.
    scale_y (float): Scaling factor for Y-position.
    font_name (str): Font name for the text.
    font_interpreter (str): Interpreter type (not applicable in matplotlib, kept for compatibility).
    font_size (int): Font size for the text.
    """
    import logging
    # Suppress just the 'findfont' warnings from matplotlib's font manager
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

    if (version_str):

        x_lim = ax.get_xlim()
        y_lim = ax.get_ylim()

        eps = 1.e-10
        if plot_type == 'loglog':
            x_lim = np.maximum(eps,x_lim)
            y_lim = np.maximum(eps,y_lim)
            x_pos = 10**(np.log10(x_lim[0]) + scale_x * (np.log10(x_lim[1]) - np.log10(x_lim[0])))
            y_pos = 10**(np.log10(y_lim[0]) + scale_y * (np.log10(y_lim[1]) - np.log10(y_lim[0])))
        elif plot_type == 'semilogx':
            x_lim = np.maximum(eps,x_lim)
            x_pos = 10**(np.log10(x_lim[0]) + scale_x * (np.log10(x_lim[1]) - np.log10(x_lim[0])))
            y_pos = y_lim[0] + scale_y * (y_lim[1] - y_lim[0])
        elif plot_type == 'semilogy':
            x_pos = x_lim[0] + scale_x * (x_lim[1] - x_lim[0])
            y_lim = np.maximum(eps,y_lim)
            y_pos = 10**(np.log10(y_lim[0]) + scale_y * (np.log10(y_lim[1]) - np.log10(y_lim[0])))
        else:
            x_pos = x_lim[0] + scale_x * (x_lim[1] - x_lim[0])
            y_pos = y_lim[0] + scale_y * (y_lim[1] - y_lim[0])

        ax.text(x_pos, y_pos, version_str, fontsize=font_size, fontname=font_name, verticalalignment='bottom', horizontalalignment='right')



def get_plot_style(style="fds"):
    """
    Returns a dictionary of plot style parameters based on the specified style.

    Parameters:
    - style (str): The style to use ('fds', 'paper', etc.). Default is 'fds'.

    Returns:
    - dict: A dictionary containing plot style parameters.
    """
    import logging
    # Suppress just the 'findfont' warnings from matplotlib's font manager
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

    if style == "fds":
        return {
            # Font properties
            "Font_Name": "Times",
            "Font_Interpreter": "TeX",
            "Key_Font_Size": 12,
            "Title_Font_Size": 16,
            "Label_Font_Size": 16,
            "Scat_Title_Font_Size": 14,
            "Scat_Label_Font_Size": 14,
            "Marker_Size": 4,
            "D1_Marker_Size": 4,
            "D2_Marker_Size": 4,

            # Line properties
            "Line_Width": 1.0,
            "D1_Line_Width": 1.0,
            "D2_Line_Width": 1.0,

            # Plot properties
            "Plot_Units": "inches",
            "Plot_Width": 5.0,
            "Plot_Height": 3.4,
            "Plot_X": 1.2,
            "Plot_Y": 0.8,
            "Scat_Plot_Width": 4.75,
            "Scat_Plot_Height": 4.75,
            "Scat_Plot_X": 1.00,
            "Scat_Plot_Y": 0.75,
            "Subtitle_Text_Offset": 0.05,
            "VerStr_Scale_X": 0.60,
            "VerStr_Scale_Y": 1.05,

            # Paper properties
            "Paper_Units": "inches",
            "Paper_Width": 6.5,
            "Paper_Height": 4.5,
            "Scat_Paper_Height": 6.0,
            "Scat_Paper_Width": 6.0,

            # Print properties
            "Figure_Visibility": "off",
            "Image_File_Type": "-dpdf",
        }

    elif style == "paper":
        return {
            # Font properties
            "Font_Name": "Helvetica",
            "Font_Interpreter": "TeX",
            "Key_Font_Size": 16,
            "Title_Font_Size": 20,
            "Label_Font_Size": 20,
            "Scat_Title_Font_Size": 14,
            "Scat_Label_Font_Size": 14,
            "Marker_Size": 10,
            "D1_Marker_Size": 10,
            "D2_Marker_Size": 10,

            # Line properties
            "Line_Width": 1.0,
            "D1_Line_Width": 1.0,
            "D2_Line_Width": 1.0,

            # Plot properties
            "Plot_Units": "normalized",
            "Plot_X": 0.1500,
            "Plot_Y": 0.1500,
            "Plot_Width": 0.7750,
            "Plot_Height": 0.8150 * 0.95,  # Adjusted for exponential y-axis tick labels
            "Scat_Plot_Width": 4.75,
            "Scat_Plot_Height": 4.75,
            "Scat_Plot_X": 0.75,
            "Scat_Plot_Y": 0.75,
            "Subtitle_Text_Offset": 0.05,
            "VerStr_Scale_X": 0.60,
            "VerStr_Scale_Y": 1.05,

            # Paper properties
            "Paper_Units": "inches",
            "Paper_Width": 8.0,
            "Paper_Height": 6.0,
            "Scat_Paper_Height": 6.0,
            "Scat_Paper_Width": 6.0,

            # Print properties
            "Figure_Visibility": "on",
            "Image_File_Type": "-dpdf",
        }

    else:
        raise ValueError(f"Unknown style '{style}'. Please choose 'fds' or 'paper'.")


def define_plot_parameters(C, irow):
    import numpy as np
    import re

    class plot_parameters:
        def __init__(self):
            self.switch_id            = C.values[irow,C.columns.get_loc('switch_id')]
            self.Dataname             = C.values[irow,C.columns.get_loc('Dataname')]
            self.VerStr_Filename      = C.values[irow,C.columns.get_loc('VerStr_Filename')]
            self.d1_Filename          = C.values[irow,C.columns.get_loc('d1_Filename')]
            self.d1_Col_Name_Row      = C.values[irow,C.columns.get_loc('d1_Col_Name_Row')]
            self.d1_Data_Row          = C.values[irow,C.columns.get_loc('d1_Data_Row')]
            self.d1_Ind_Col_Name      = C.values[irow,C.columns.get_loc('d1_Ind_Col_Name')]
            self.d1_Dep_Col_Name      = C.values[irow,C.columns.get_loc('d1_Dep_Col_Name')]
            self.d1_Key               = C.values[irow,C.columns.get_loc('d1_Key')]
            self.d1_Style             = C.values[irow,C.columns.get_loc('d1_Style')]
            self.d1_Start             = C.values[irow,C.columns.get_loc('d1_Start')]
            self.d1_End               = C.values[irow,C.columns.get_loc('d1_End')]
            self.d1_Tick              = C.values[irow,C.columns.get_loc('d1_Tick')]
            self.d1_Comp_Start        = C.values[irow,C.columns.get_loc('d1_Comp_Start')]
            self.d1_Comp_End          = C.values[irow,C.columns.get_loc('d1_Comp_End')]
            self.d1_Dep_Comp_Start    = C.values[irow,C.columns.get_loc('d1_Dep_Comp_Start')]
            self.d1_Dep_Comp_End      = C.values[irow,C.columns.get_loc('d1_Dep_Comp_End')]
            self.d1_Initial_Value     = C.values[irow,C.columns.get_loc('d1_Initial_Value')]
            self.d2_Filename          = C.values[irow,C.columns.get_loc('d2_Filename')]
            self.d2_Col_Name_Row      = C.values[irow,C.columns.get_loc('d2_Col_Name_Row')]
            self.d2_Data_Row          = C.values[irow,C.columns.get_loc('d2_Data_Row')]
            self.d2_Ind_Col_Name      = C.values[irow,C.columns.get_loc('d2_Ind_Col_Name')]
            self.d2_Dep_Col_Name      = C.values[irow,C.columns.get_loc('d2_Dep_Col_Name')]
            self.d2_Key               = C.values[irow,C.columns.get_loc('d2_Key')]
            self.d2_Style             = C.values[irow,C.columns.get_loc('d2_Style')]
            self.d2_Start             = C.values[irow,C.columns.get_loc('d2_Start')]
            self.d2_End               = C.values[irow,C.columns.get_loc('d2_End')]
            self.d2_Tick              = C.values[irow,C.columns.get_loc('d2_Tick')]
            self.d2_Comp_Start        = C.values[irow,C.columns.get_loc('d2_Comp_Start')]
            self.d2_Comp_End          = C.values[irow,C.columns.get_loc('d2_Comp_End')]
            self.d2_Dep_Comp_Start    = C.values[irow,C.columns.get_loc('d2_Dep_Comp_Start')]
            self.d2_Dep_Comp_End      = C.values[irow,C.columns.get_loc('d2_Dep_Comp_End')]
            self.d2_Initial_Value     = C.values[irow,C.columns.get_loc('d2_Initial_Value')]
            self.Plot_Title           = C.values[irow,C.columns.get_loc('Plot_Title')]
            self.Ind_Title            = C.values[irow,C.columns.get_loc('Ind_Title')]
            self.Dep_Title            = C.values[irow,C.columns.get_loc('Dep_Title')]
            self.Min_Ind              = C.values[irow,C.columns.get_loc('Min_Ind')]
            self.Max_Ind              = C.values[irow,C.columns.get_loc('Max_Ind')]
            self.Scale_Ind            = C.values[irow,C.columns.get_loc('Scale_Ind')]
            self.Min_Dep              = C.values[irow,C.columns.get_loc('Min_Dep')]
            self.Max_Dep              = C.values[irow,C.columns.get_loc('Max_Dep')]
            self.Scale_Dep            = C.values[irow,C.columns.get_loc('Scale_Dep')]
            self.Flip_Axis            = C.values[irow,C.columns.get_loc('Flip_Axis')]
            self.Title_Position       = C.values[irow,C.columns.get_loc('Title_Position')]
            self.Key_Position         = C.values[irow,C.columns.get_loc('Key_Position')]
            self.Legend_XYWidthHeight = C.values[irow,C.columns.get_loc('Legend_XYWidthHeight')]
            self.Paper_Width_Factor   = C.values[irow,C.columns.get_loc('Paper_Width_Factor')]
            self.Plot_Type            = C.values[irow,C.columns.get_loc('Plot_Type')]
            self.Plot_Filename        = C.values[irow,C.columns.get_loc('Plot_Filename')]
            self.Quantity             = C.values[irow,C.columns.get_loc('Quantity')]
            self.Metric               = C.values[irow,C.columns.get_loc('Metric')]
            self.Error_Tolerance      = C.values[irow,C.columns.get_loc('Error_Tolerance')]
            self.Group_Key_Label      = C.values[irow,C.columns.get_loc('Group_Key_Label')]
            self.Group_Style          = C.values[irow,C.columns.get_loc('Group_Style')]
            self.Fill_Color           = C.values[irow,C.columns.get_loc('Fill_Color')]
            self.Font_Interpreter     = C.values[irow,C.columns.get_loc('Font_Interpreter')]

        def __repr__(self):
            return str(self.__dict__)

    inst = plot_parameters()

    specials = {
        "&": r"\&", "%": r"\%", "_": r"\_", "#": r"\#",
        "$": r"\$", "{": r"\{", "}": r"\}", "^": r"\^{}", "~": r"\~{}",
    }

    def sanitize(text: str) -> str:
        parts = re.split(r"(\$.*?\$)", text)
        sanitized = []
        for part in parts:
            if part.startswith("$") and part.endswith("$"):
                sanitized.append(part)  # math untouched
            else:
                s = part
                for k, v in specials.items():
                    s = s.replace(k, v)
                sanitized.append(s)
        return "".join(sanitized)

    def safe_strip(val):
        if isinstance(val, str):
            return val.strip()
        return ""  # or return val if you prefer to keep None

    # Explicit sanitization of only the human-facing fields
    inst.Plot_Title      = sanitize(safe_strip(inst.Plot_Title))
    inst.Ind_Title       = sanitize(safe_strip(inst.Ind_Title))
    inst.Dep_Title       = sanitize(safe_strip(inst.Dep_Title))
    inst.Quantity        = sanitize(safe_strip(inst.Quantity))
    inst.Metric          = sanitize(safe_strip(inst.Metric))
    inst.Group_Key_Label = sanitize(safe_strip(inst.Group_Key_Label))
    inst.d1_Key          = sanitize(safe_strip(inst.d1_Key))
    inst.d2_Key          = sanitize(safe_strip(inst.d2_Key))

    return inst


def matlab_legend_to_matplotlib(position):
    """
    Convert a MATLAB legend position string to the corresponding matplotlib location string.

    Parameters:
        position (str): MATLAB-style legend position (e.g., 'northeast', 'southwestoutside')

    Returns:
        str: Matplotlib-compatible legend location (e.g., 'upper right')
    """
    mapping = {
        'north': 'upper center',
        'south': 'lower center',
        'east': 'center right',
        'west': 'center left',
        'northeast': 'upper right',
        'southeast': 'lower right',
        'southwest': 'lower left',
        'northwest': 'upper left',
        'northeastoutside': 'center left',    # rough equivalent
        'northwestoutside': 'center right',
        'southeastoutside': 'center left',
        'southwestoutside': 'center right',
        'best': 'best'
    }

    if not isinstance(position, str):
        return 'best'

    return mapping.get(position.strip().lower(), 'best')


def define_qrow_variables(Q, j):
    """
    Define scatterplot parameters from the Scatterplot_Inputs.csv row j.

    Mirrors the MATLAB 'define_qrow_variables.m' behavior and returns
    a simple Python object with all scatterplot configuration fields.

    Parameters
    ----------
    Q : pandas.DataFrame
        The scatterplot input file loaded by pandas.read_csv().
    j : int
        Row index (0-based).

    Returns
    -------
    q : object
        Object with attributes corresponding to scatterplot input fields.
    """
    import re

    class qrow:
        def __init__(self):
            self.Scatter_Plot_Title  = Q.loc[j, "Scatter_Plot_Title"]
            self.Ind_Title           = Q.loc[j, "Ind_Title"]
            self.Dep_Title           = Q.loc[j, "Dep_Title"]
            self.Plot_Min            = Q.loc[j, "Plot_Min"]
            self.Plot_Max            = Q.loc[j, "Plot_Max"]
            self.Title_Position      = Q.loc[j, "Title_Position"]
            self.Key_Position        = Q.loc[j, "Key_Position"]
            self.Paper_Width_Factor  = Q.loc[j, "Paper_Width_Factor"]
            self.Sigma_E             = Q.loc[j, "Sigma_E"]
            self.Weight_Data         = Q.loc[j, "Weight_Data"]
            self.Plot_Type           = Q.loc[j, "Plot_Type"]
            self.Plot_Filename       = Q.loc[j, "Plot_Filename"]

        def __repr__(self):
            return str(self.__dict__)

    q = qrow()

    # --- Sanitization helpers (same style as define_plot_parameters) ---
    specials = {
        "&": r"\&", "%": r"\%", "_": r"\_", "#": r"\#",
        "$": r"\$", "{": r"\{", "}": r"\}", "^": r"\^{}", "~": r"\~{}",
    }

    def sanitize(text: str) -> str:
        if not isinstance(text, str):
            return text
        parts = re.split(r"(\$.*?\$)", text)
        sanitized = []
        for part in parts:
            if part.startswith("$") and part.endswith("$"):
                sanitized.append(part)
            else:
                s = part
                for k, v in specials.items():
                    s = s.replace(k, v)
                sanitized.append(s)
        return "".join(sanitized)

    def safe_strip(val):
        if isinstance(val, str):
            return val.strip()
        return val

    # Sanitize only human-readable fields
    q.Scatter_Plot_Title = sanitize(safe_strip(q.Scatter_Plot_Title))
    q.Ind_Title = sanitize(safe_strip(q.Ind_Title))
    q.Dep_Title = sanitize(safe_strip(q.Dep_Title))
    q.Plot_Filename = safe_strip(q.Plot_Filename)
    q.Key_Position = safe_strip(q.Key_Position)
    q.Plot_Type = safe_strip(q.Plot_Type)

    # Parse numeric fields
    def to_float(val):
        try:
            return float(val)
        except Exception:
            return np.nan

    q.Plot_Min = to_float(q.Plot_Min)
    q.Plot_Max = to_float(q.Plot_Max)
    q.Paper_Width_Factor = to_float(q.Paper_Width_Factor)
    q.Sigma_E = to_float(q.Sigma_E)

    # Parse Title_Position as [x, y] floats
    if isinstance(q.Title_Position, str):
        try:
            vals = [float(x) for x in q.Title_Position.split()]
            if len(vals) == 2:
                q.Title_Position = vals
            else:
                q.Title_Position = [0.03, 0.95]
        except Exception:
            q.Title_Position = [0.03, 0.95]
    else:
        q.Title_Position = [0.03, 0.95]

    # Normalize Weight_Data (yes/no → bool)
    if isinstance(q.Weight_Data, str):
        q.Weight_Data = q.Weight_Data.strip().lower() == "yes"
    else:
        q.Weight_Data = bool(q.Weight_Data)

    return q


def scatplot(saved_data, drange, **kwargs):
    import os
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    # --- Inputs ---
    manuals_dir = kwargs.get('Manuals_Dir', '')
    scatterplot_inputs_file = kwargs.get('Scatterplot_Inputs_File', '')
    stats_output = kwargs.get('Stats_Output', 'Validation')
    scatterplot_dir = kwargs.get('Scatterplot_Dir', '')
    verbose = kwargs.get('verbose', True)

    # --- Prepare paths ---
    scatterplot_path = scatterplot_inputs_file  # use exactly what the caller passed
    if not os.path.exists(scatterplot_path):
        raise FileNotFoundError(f"[scatplot] Cannot find Scatterplot_Inputs_File: {scatterplot_path}")
    
    df = pd.read_csv(scatterplot_path, sep=',', engine='python', quotechar='"')

    if verbose:
        print(f"[scatplot] Loaded {len(df)} lines from scatterplot inputs")

    # --- Unpack dataplot output ---
    (Save_Quantity,
     Save_Group_Style,
     Save_Fill_Color,
     Save_Group_Key_Label,
     Save_Measured_Metric,
     Save_Predicted_Metric,
     Save_Dataname,
     Save_Plot_Filename,
     Save_Dep_Title,
     Save_Error_Tolerance,
     Save_Metric_Type,
     Save_Measured_Quantity,
     Save_Predicted_Quantity) = saved_data

    # Ensure output directory exists
    os.makedirs(scatterplot_dir, exist_ok=True)

    # --- Loop over scatterplot CSV rows ---
    for idx, row in df.iterrows():
        title = row['Scatter_Plot_Title']
        ind_title = row['Ind_Title']
        dep_title = row['Dep_Title']
        plot_min = float(row['Plot_Min'])
        plot_max = float(row['Plot_Max'])
        plot_type = row['Plot_Type'].strip().lower()
        plot_filename = os.path.join(manuals_dir, row['Plot_Filename'] + '.pdf')

        if verbose:
            print(f"[scatplot] Generating scatter plot: {title}")

        # --- Find dataplot rows matching this title ---
        match_indices = [i for i, q in enumerate(Save_Quantity)
                         if title.split(',')[0].strip().lower() in str(q).lower()]

        if not match_indices:
            print(f"[scatplot] No dataplot entries matched '{title}' — skipping.")
            continue

        # --- Collect measured/predicted data ---
        measured_vals = []
        predicted_vals = []
        for i in match_indices:
            m = np.atleast_1d(Save_Measured_Metric[i]).astype(float)
            p = np.atleast_1d(Save_Predicted_Metric[i]).astype(float)
            # Filter out invalids
            mask = np.isfinite(m) & np.isfinite(p)
            measured_vals.extend(m[mask])
            predicted_vals.extend(p[mask])

        measured_vals = np.array(measured_vals)
        predicted_vals = np.array(predicted_vals)

        if len(measured_vals) == 0:
            print(f"[scatplot] Warning: No valid data points for '{title}'")
            continue

        # --- Create scatter plot ---
        plt.figure(figsize=(5,5))
        if plot_type == 'loglog':
            plt.loglog(measured_vals, predicted_vals, 'o', markersize=5, alpha=0.6)
        else:
            plt.plot(measured_vals, predicted_vals, 'o', markersize=5, alpha=0.6)

        # 1:1 line
        plt.plot([plot_min, plot_max], [plot_min, plot_max], 'k--', linewidth=1)

        plt.xlim([plot_min, plot_max])
        plt.ylim([plot_min, plot_max])
        plt.xlabel(ind_title)
        plt.ylabel(dep_title)
        plt.title(title)
        plt.grid(True, which='both', linestyle=':')

        plt.tight_layout()
        os.makedirs(os.path.dirname(plot_filename), exist_ok=True)
        plt.savefig(plot_filename)
        plt.close()


    # --- Placeholder output_stats table ---
    # This will later be replaced with computed statistics per quantity
    # For now, collect each scatterplot title with basic counts
    output_stats = []
    output_stats.append([
        "Quantity", "Datasets", "Points", "Sigma_E", "Sigma_M", "Bias"
    ])

    for idx, row in df.iterrows():
        title = row['Scatter_Plot_Title']
        match_indices = [i for i, q in enumerate(Save_Quantity)
                         if title.split(',')[0].strip().lower() in str(q).lower()]
        if not match_indices:
            continue
        num_datasets = len(match_indices)
        num_points = sum(
            len(np.atleast_1d(Save_Measured_Metric[i])) for i in match_indices
        )
        sigma_e = np.nanstd([Save_Measured_Metric[i] for i in match_indices])
        sigma_m = np.nanstd([Save_Predicted_Metric[i] for i in match_indices])
        bias = np.nanmean([
            np.nanmean(Save_Predicted_Metric[i]) - np.nanmean(Save_Measured_Metric[i])
            for i in match_indices
        ])
        output_stats.append([
            title, num_datasets, num_points, sigma_e, sigma_m, bias
        ])


    # --- Stats output placeholder ---
    stats_csv = os.path.join(scatterplot_dir, f"ScatterPlot_Stats_{stats_output}.csv")
    tex_out   = os.path.join(scatterplot_dir, f"ScatterPlot_Tables_{stats_output}.tex")
    hist_out  = os.path.join(scatterplot_dir, f"ScatterPlot_Histograms_{stats_output}.tex")

    statistics_output(
        Stats_Output=stats_output,
        output_stats=output_stats,
        Output_File=stats_csv,
        Statistics_Tex_Output=tex_out,
        Histogram_Tex_Output=hist_out,
        Output_Histograms=[]
    )


    print("[scatplot] Completed successfully.")



def statistics_output(
    Stats_Output,
    output_stats,
    Output_File,
    Statistics_Tex_Output=None,
    Histogram_Tex_Output=None,
    Output_Histograms=None,
):
    """
    Python translation of statistics_output.m (K. Overholt, 2013)
    Modernized using pandas for CSV I/O.

    Parameters
    ----------
    Stats_Output : str
        'Validation', 'Verification', or 'None'
    output_stats : list of lists, numpy array, or pandas DataFrame
        Statistical results to write
    Output_File : str
        Path to write CSV file
    Statistics_Tex_Output : str, optional
        Path to write LaTeX table file
    Histogram_Tex_Output : str, optional
        Path to write LaTeX histogram file
    Output_Histograms : list, optional
        List of histogram figure filenames
    """
    import os
    import numpy as np
    import pandas as pd

    if Stats_Output == "None":
        print("[statistics_output] Skipping (Stats_Output=None)")
        return

    # --- normalize to DataFrame ---
    if not isinstance(output_stats, pd.DataFrame):
        try:
            output_stats = pd.DataFrame(output_stats)
        except Exception as e:
            print(f"[statistics_output] Could not convert output_stats to DataFrame: {e}")
            return

    if output_stats.empty:
        print("[statistics_output] No stats to write")
        return

    # --- Write CSV ---
    os.makedirs(os.path.dirname(Output_File), exist_ok=True)
    output_stats.to_csv(Output_File, index=False, quoting=1)  # quote all non-numerics
    print(f"[statistics_output] Wrote CSV: {Output_File}")

    # --- Write LaTeX for Verification ---
    if Stats_Output == "Verification" and Statistics_Tex_Output:
        with open(Statistics_Tex_Output, "w") as fid:
            fid.write("\\scriptsize\n")
            fid.write("\\begin{longtable}{|p{2.5in}|l|p{1in}|l|p{1in}|l|l|l|l|l|}\n")
            fid.write("\\hline\n")
            fid.write(
                "Case Name & Section & Expected & Expected & Predicted & Predicted & "
                "Type of & Error & Error & Within \\\\\n"
            )
            fid.write(
                "          &         & Quantity & Value & Quantity & Value & Error &  "
                "& Tolerance & Tol. \\\\ \\hline \\hline\n"
            )
            fid.write("\\endfirsthead\n\\hline\n")
            fid.write(
                "Case Name & Section & Expected & Expected & Predicted & Predicted & "
                "Type of & Error & Error & Within \\\\\n"
            )
            fid.write(
                "          &         & Quantity & Value & Quantity & Value & Error &  "
                "& Tolerance & Tol. \\\\ \\hline \\hline\n"
            )
            fid.write("\\endhead\n\\hline\n\\endfoot\n\\hline\n\\endlastfoot\n")

            m = output_stats.iloc[1:].sort_values(by=output_stats.columns[0])
            for _, r in m.iterrows():
                try:
                    case = f"\\lstinline[basicstyle=\\scriptsize\\ttfamily]!{r[2]}!"
                    section = f"\\ref{{{r[2]}}}"
                    exp_q = f"\\lstinline[basicstyle=\\scriptsize\\ttfamily]!{r[4]}!"
                    pred_q = f"\\lstinline[basicstyle=\\scriptsize\\ttfamily]!{r[6]}!"
                    fid.write(
                        f"{case} & {section} & {exp_q} & {float(r[5]):1.2e} & "
                        f"{pred_q} & {float(r[7]):1.2e} & "
                        f"{str(r[9]).replace(' Error','')} & "
                        f"{float(r[10]):1.2e} & {float(r[11]):1.2e} & {r[12]} \\\\\n"
                    )
                except Exception as e:
                    print(f"[statistics_output] Skipped Verification row due to error: {e}")

            fid.write("\\end{longtable}\n\\normalsize\n")
        print(f"[statistics_output] Wrote LaTeX Verification table: {Statistics_Tex_Output}")

    # --- Write LaTeX for Validation ---
    if Stats_Output == "Validation" and Statistics_Tex_Output:
        with open(Statistics_Tex_Output, "w") as fid:
            fid.write("\\begin{longtable}[c]{|l|c|c|c|c|c|c|}\n")
            fid.write(
                "\\caption[Summary statistics]{Summary statistics for all quantities of interest}\n"
            )
            fid.write("\\label{summary_stats}\n\\\\ \\hline\n")
            fid.write(
                "Quantity & Section & Datasets & Points & "
                "$\\widetilde{\\sigma}_{\\rm E}$ & $\\widetilde{\\sigma}_{\\rm M}$ & Bias \\\\ \\hline \\hline\n"
            )
            fid.write("\\endfirsthead\n\\hline\n")
            fid.write(
                "Quantity & Section & Datasets & Points & "
                "$\\widetilde{\\sigma}_{\\rm E}$ & $\\widetilde{\\sigma}_{\\rm M}$ & Bias \\\\ \\hline \\hline\n"
            )
            fid.write("\\endhead\n")

            for _, r in output_stats.iloc[1:].iterrows():
                try:
                    sigma_e = float(r[3])
                    if sigma_e < 0:
                        continue
                    quantity = r[0]
                    section = f"\\ref{{{quantity}}}"
                    fid.write(
                        f"{quantity} & {section} & {int(r[1])} & {int(r[2])} & "
                        f"{sigma_e:0.2f} & {float(r[4]):0.2f} & {float(r[5]):0.2f} \\\\ \\hline\n"
                    )
                except Exception as e:
                    print(f"[statistics_output] Skipped Validation row due to error: {e}")

            fid.write("\\end{longtable}\n")
        print(f"[statistics_output] Wrote LaTeX Validation table: {Statistics_Tex_Output}")

    # --- Histogram LaTeX (optional) ---
    if Stats_Output == "Validation" and Output_Histograms:
        os.makedirs(os.path.dirname(Histogram_Tex_Output), exist_ok=True)
        with open(Histogram_Tex_Output, "w") as fid:
            num_hist = len(Output_Histograms)
            pages = int(np.ceil(num_hist / 8))
            for i in range(pages):
                fid.write("\\begin{figure}[p]\n")
                fid.write("\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}r}\n")
                for j in range(i * 8, min((i + 1) * 8, num_hist)):
                    line_ending = "&" if j % 2 == 0 else "\\\\"
                    fid.write(
                        f"\\includegraphics[height=2.2in]"
                        f"{{SCRIPT_FIGURES/ScatterPlots/{Output_Histograms[j]}}} {line_ending}\n"
                    )
                fid.write("\\end{tabular*}\n")
                fid.write(f"\\label{{Histogram_{i+1}}}\n\\end{figure}\n\n")
        print(f"[statistics_output] Wrote LaTeX histograms: {Histogram_Tex_Output}")
