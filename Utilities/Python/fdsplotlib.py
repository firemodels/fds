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
import numpy as np
import pandas as pd


def dataplot(config_filename,**kwargs):

    import os
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib import rc

    import logging
    # Suppress just the 'findfont' warnings from matplotlib's font manager
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

    # defaults
    revision  = ''
    configdir = ''
    expdir = ''
    cmpdir = ''
    pltdir = ''
    close_figs = False
    verbose = False
    plot_list = ['all']
    plot_range = range(10000)

    if kwargs.get('revision'):
        revision = kwargs.get('revision')

    if kwargs.get('expdir'):
        expdir = kwargs.get('expdir')

    if kwargs.get('cmpdir'):
        cmpdir = kwargs.get('cmpdir')

    if kwargs.get('pltdir'):
        pltdir = kwargs.get('pltdir')

    if kwargs.get('close_figs'):
        close_figs = kwargs.get('close_figs')

    if kwargs.get('verbose'):
        verbose = kwargs.get('verbose')

    if kwargs.get('plot_list'):
        plot_list = kwargs.get('plot_list')

    if kwargs.get('plot_range'):
        plot_range_in = kwargs.get('plot_range')
        plot_range = range(plot_range_in[0]-2,plot_range_in[-1]-1)

    # read the config file
    df = pd.read_csv(configdir+config_filename, sep=',', engine='python', comment='#', quotechar='"')
    C = df.where(pd.notnull(df), None)

    Plot_Filename_Last = None
    d1_Key_Last = None
    f_Last = plt.figure()

    # loop over the rows of the config file
    for irow in C.index:

        if irow not in plot_range:
            continue

        # try:

        # define plot parameters and return them in an object called pp
        pp = define_plot_parameters(C,irow)

        # print(pp.__dict__) # helpful for debug

        if 'all' not in plot_list:
            if pp.Dataname not in plot_list:
                continue

        if pp.Plot_Filename!=Plot_Filename_Last:

            if verbose:
                print('Generating plot ' + pltdir + pp.Plot_Filename + '...')

            if close_figs:
                plt.close('all')

            # read data from exp file
            # set header to the row where column names are stored (Python is 0 based)
            E = pd.read_csv(expdir+pp.d1_Filename, header=int(pp.d1_Col_Name_Row-1), sep=',', engine='python', comment='#', quotechar='"')

            x = E[pp.d1_Ind_Col_Name].values[:].astype(float)
            # y = E[pp.d1_Dep_Col_Name].values[:].astype(float)
            col_names = [c.strip() for c in pp.d1_Dep_Col_Name.split('|')]
            # print(col_names)
            y = E[col_names].values.astype(float)

            for i, label in enumerate(col_names):
                # plot the exp data
                f = plot_to_fig(x_data=x, y_data=y[:, i],
                    data_label=pp.d1_Key,
                    x_label=pp.Ind_Title,
                    y_label=pp.Dep_Title,
                    marker_style=pp.d1_Style,
                    x_min=pp.Min_Ind,x_max=pp.Max_Ind,
                    y_min=pp.Min_Dep,y_max=pp.Max_Dep,
                    legend_location=pp.Key_Position
                    )

            # plt.figure(f.number) # make figure current
            # plt.show()
        else:
            f = f_Last

            if pp.d1_Key!=d1_Key_Last:

                # read data from exp file
                # set header to the row where column names are stored (Python is 0 based)
                E = pd.read_csv(expdir+pp.d1_Filename, header=int(pp.d1_Col_Name_Row-1), sep=',', engine='python', comment='#', quotechar='"')
                x = E[pp.d1_Ind_Col_Name].values[:].astype(float)
                y = E[pp.d1_Dep_Col_Name].values[:].astype(float)

                # plot the exp data
                f = plot_to_fig(x_data=x, y_data=y,
                    figure_handle=f,
                    data_label=pp.d1_Key,
                    x_label=pp.Ind_Title,
                    y_label=pp.Dep_Title,
                    marker_style=pp.d1_Style,
                    x_min=pp.Min_Ind,x_max=pp.Max_Ind,
                    y_min=pp.Min_Dep,y_max=pp.Max_Dep,
                    legend_location=pp.Key_Position
                    )

        # get the model results
        M = pd.read_csv(cmpdir+pp.d2_Filename, header=int(pp.d2_Col_Name_Row-1), sep=',', engine='python', comment='#', quotechar='"')
        x = M[pp.d2_Ind_Col_Name].values[:].astype(float)
        # y = M[pp.d2_Dep_Col_Name].values[:].astype(float)
        col_names = [c.strip() for c in pp.d2_Dep_Col_Name.split('|')]
        y = M[col_names].values.astype(float)

        version_string = revision
        if (pp.VerStr_Filename):
            file1 = open(cmpdir+pp.VerStr_Filename,"r")
            Lines = file1.readlines()
            version_string = Lines[0].strip()
            file1.close()

        for i, label in enumerate(col_names):
            f = plot_to_fig(x_data=x, y_data=y[:, i],
                revision_label=version_string,
                figure_handle=f,
                x_label=pp.Ind_Title,
                y_label=pp.Dep_Title,
                data_label=pp.d2_Key,
                line_style=pp.d2_Style,
                x_min=pp.Min_Ind,x_max=pp.Max_Ind,
                y_min=pp.Min_Dep,y_max=pp.Max_Dep,
                legend_location=pp.Key_Position,
                plot_title=pp.Plot_Title
                )

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
    default_ticklabel_fontsize = plot_style["Label_Font_Size"]
    default_axeslabel_fontsize = plot_style["Label_Font_Size"]
    default_legend_fontsize = plot_style["Key_Font_Size"]
    default_title_fontsize = plot_style["Title_Font_Size"]
    default_markevery = 1
    markerfacecolor = 'none'
    markeredgecolor = 'black'
    markeredgewidth = 1
    marker = None
    markersize = 5
    linestyle = '-'
    linewidth = 1
    color = 'black'
    ###############################

    if kwargs.get('figure_size'):
        figure_size=kwargs.get('figure_size')
    else:
        figure_size=default_figure_size

    # if figure handle is passed, append to current figure, else generate a new figure
    if kwargs.get('figure_handle'):
        fig = kwargs.get('figure_handle')
        ax = fig.axes[0]
        plt.figure(fig.number)
    else:
        fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, gridspec_kw={'hspace': 0, 'wspace': 0}, figsize=figure_size)

    # select plot type
    if kwargs.get('plot_type'):
        plot_type=kwargs.get('plot_type')
    else:
        plot_type='linear'

    # convert matlab styles to matplotlib
    if kwargs.get('marker_style'):
        style = kwargs.get('marker_style')

        color,marker,linestyle = parse_matlab_style(style)


    # other plot parameters
    if kwargs.get('data_markevery'):
        markevery = kwargs.get('data_markevery')
    else:
        markevery = default_markevery


    # generate the main x,y plot
    if plot_type=='linear':
        ax.plot(x_data,y_data,
            markevery=markevery,
            label=kwargs.get('data_label'),
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
            label=kwargs.get('data_label'),
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
            label=kwargs.get('data_label'),
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
            label=kwargs.get('data_label'),
            markerfacecolor=markerfacecolor,
            markeredgecolor=color,
            markeredgewidth=markeredgewidth,
            marker=marker,
            markersize=markersize,
            linestyle=linestyle,
            linewidth=linewidth,
            color=color)

    # if error range is passed, add it to the plot
    if kwargs.get('y_error_absolute') and not kwargs.get('y_error_relative'):
        if kwargs.get('y_error_absolute')>0.:
            ax.fill_between(x_data,y_data-kwargs.get('y_error_absolute'),y_data+kwargs.get('y_error_absolute'),
                alpha=0.1,color=kwargs.get('marker_edge_color'))

    if kwargs.get('y_error_relative') and not kwargs.get('y_error_absolute'):
        if kwargs.get('y_error_relative')>0.:
            ax.fill_between(x_data,y_data*(1.-kwargs.get('y_error_relative')),y_data*(1.+kwargs.get('y_error_relative')),
                alpha=0.1,color=kwargs.get('marker_edge_color'))

    if kwargs.get('y_error_relative') and kwargs.get('y_error_absolute'):
        if kwargs.get('y_error_relative')>0.:
            ax.fill_between(x_data,y_data*(1.-kwargs.get('y_error_relative'))-kwargs.get('y_error_absolute'),y_data*(1.+kwargs.get('y_error_relative'))+kwargs.get('y_error_absolute'),
                alpha=0.1,color=kwargs.get('marker_edge_color'))

    try:
        y_error = kwargs.get('y_error_vector')
        if len(y_data)==len(y_error):
            ax.fill_between(x_data,y_data-y_error,y_data+y_error,
                alpha=0.1,color=kwargs.get('marker_edge_color'))
    except:
        y_error = 0.

    if kwargs.get('ticklabel_fontsize'):
        ticklabel_fontsize=kwargs.get('ticklabel_fontsize')
    else:
        ticklabel_fontsize=default_ticklabel_fontsize

    plt.setp( ax.xaxis.get_majorticklabels(), rotation=0, fontsize=ticklabel_fontsize )
    plt.setp( ax.yaxis.get_majorticklabels(), rotation=0, fontsize=ticklabel_fontsize )

    if kwargs.get('axeslabel_fontsize'):
        axeslabel_fontsize=kwargs.get('axeslabel_fontsize')
    else:
        axeslabel_fontsize=default_axeslabel_fontsize

    plt.xlabel(kwargs.get('x_label'), fontsize=axeslabel_fontsize)
    plt.ylabel(kwargs.get('y_label'), fontsize=axeslabel_fontsize)

    if kwargs.get('legend_fontsize'):
        legend_fontsize=kwargs.get('legend_fontsize')
    else:
        legend_fontsize=default_legend_fontsize

    if kwargs.get('legend_location')=='outside':
        plt.legend(fontsize=legend_fontsize,bbox_to_anchor=(1,1),loc='upper left',framealpha=kwargs.get('legend_framealpha'))
    else:
        if kwargs.get('show_legend'):
            plt.legend(fontsize=legend_fontsize,loc=kwargs.get('legend_location'),framealpha=kwargs.get('legend_framealpha'))

    # plot titles
    if kwargs.get('title_fontsize'):
        title_fontsize=kwargs.get('title_fontsize')
    else:
        title_fontsize=default_title_fontsize

    # set axes and tick properties
    ymin=kwargs.get('y_min')
    ymax=kwargs.get('y_max')
    xmin=kwargs.get('x_min')
    xmax=kwargs.get('x_max')

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    if kwargs.get('revision_label'):
        add_version_string(ax, kwargs.get('revision_label'), plot_type)

    fig.tight_layout()

    return fig


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

def add_version_string(ax, version_str, plot_type='linear', scale_x=0.60, scale_y=1.02,
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

        ax.text(x_pos, y_pos, version_str, fontsize=font_size, fontname=font_name, verticalalignment='bottom')



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


def define_plot_parameters(C,irow):
    """
    gathers parameters from config file
    """
    class plot_parameters:

        def __repr__(self):
            return str(self)

        switch_id            = C.values[irow,C.columns.get_loc('switch_id')]
        Dataname             = C.values[irow,C.columns.get_loc('Dataname')]
        VerStr_Filename      = C.values[irow,C.columns.get_loc('VerStr_Filename')]
        d1_Filename          = C.values[irow,C.columns.get_loc('d1_Filename')]
        d1_Col_Name_Row      = C.values[irow,C.columns.get_loc('d1_Col_Name_Row')]
        d1_Data_Row          = C.values[irow,C.columns.get_loc('d1_Data_Row')]
        d1_Ind_Col_Name      = C.values[irow,C.columns.get_loc('d1_Ind_Col_Name')]
        d1_Dep_Col_Name      = C.values[irow,C.columns.get_loc('d1_Dep_Col_Name')]
        d1_Key               = C.values[irow,C.columns.get_loc('d1_Key')]
        d1_Style             = C.values[irow,C.columns.get_loc('d1_Style')]
        d1_Start             = C.values[irow,C.columns.get_loc('d1_Start')]
        d1_End               = C.values[irow,C.columns.get_loc('d1_End')]
        d1_Tick              = C.values[irow,C.columns.get_loc('d1_Tick')]
        d1_Comp_Start        = C.values[irow,C.columns.get_loc('d1_Comp_Start')]
        d1_Comp_End          = C.values[irow,C.columns.get_loc('d1_Comp_End')]
        d1_Dep_Comp_Start    = C.values[irow,C.columns.get_loc('d1_Dep_Comp_Start')]
        d1_Dep_Comp_End      = C.values[irow,C.columns.get_loc('d1_Dep_Comp_End')]
        d1_Initial_Value     = C.values[irow,C.columns.get_loc('d1_Initial_Value')]
        d2_Filename          = C.values[irow,C.columns.get_loc('d2_Filename')]
        d2_Col_Name_Row      = C.values[irow,C.columns.get_loc('d2_Col_Name_Row')]
        d2_Data_Row          = C.values[irow,C.columns.get_loc('d2_Data_Row')]
        d2_Ind_Col_Name      = C.values[irow,C.columns.get_loc('d2_Ind_Col_Name')]
        d2_Dep_Col_Name      = C.values[irow,C.columns.get_loc('d2_Dep_Col_Name')]
        d2_Key               = C.values[irow,C.columns.get_loc('d2_Key')]
        d2_Style             = C.values[irow,C.columns.get_loc('d2_Style')]
        d2_Start             = C.values[irow,C.columns.get_loc('d2_Start')]
        d2_End               = C.values[irow,C.columns.get_loc('d2_End')]
        d2_Tick              = C.values[irow,C.columns.get_loc('d2_Tick')]
        d2_Comp_Start        = C.values[irow,C.columns.get_loc('d2_Comp_Start')]
        d2_Comp_End          = C.values[irow,C.columns.get_loc('d2_Comp_End')]
        d2_Dep_Comp_Start    = C.values[irow,C.columns.get_loc('d2_Dep_Comp_Start')]
        d2_Dep_Comp_End      = C.values[irow,C.columns.get_loc('d2_Dep_Comp_End')]
        d2_Initial_Value     = C.values[irow,C.columns.get_loc('d2_Initial_Value')]
        Plot_Title           = C.values[irow,C.columns.get_loc('Plot_Title')]
        Ind_Title            = C.values[irow,C.columns.get_loc('Ind_Title')]
        Dep_Title            = C.values[irow,C.columns.get_loc('Dep_Title')]
        Min_Ind              = C.values[irow,C.columns.get_loc('Min_Ind')]
        Max_Ind              = C.values[irow,C.columns.get_loc('Max_Ind')]
        Scale_Ind            = C.values[irow,C.columns.get_loc('Scale_Ind')]
        Min_Dep              = C.values[irow,C.columns.get_loc('Min_Dep')]
        Max_Dep              = C.values[irow,C.columns.get_loc('Max_Dep')]
        Scale_Dep            = C.values[irow,C.columns.get_loc('Scale_Dep')]
        Flip_Axis            = C.values[irow,C.columns.get_loc('Flip_Axis')]
        Title_Position       = C.values[irow,C.columns.get_loc('Title_Position')]
        Key_Position         = C.values[irow,C.columns.get_loc('Key_Position')]
        Legend_XYWidthHeight = C.values[irow,C.columns.get_loc('Legend_XYWidthHeight')]
        Paper_Width_Factor   = C.values[irow,C.columns.get_loc('Paper_Width_Factor')]
        Plot_Type            = C.values[irow,C.columns.get_loc('Plot_Type')]
        Plot_Filename        = C.values[irow,C.columns.get_loc('Plot_Filename')]
        Quantity             = C.values[irow,C.columns.get_loc('Quantity')]
        Metric               = C.values[irow,C.columns.get_loc('Metric')]
        Error_Tolerance      = C.values[irow,C.columns.get_loc('Error_Tolerance')]
        Group_Key_Label      = C.values[irow,C.columns.get_loc('Group_Key_Label')]
        Group_Style          = C.values[irow,C.columns.get_loc('Group_Style')]
        Fill_Color           = C.values[irow,C.columns.get_loc('Fill_Color')]
        Font_Interpreter     = C.values[irow,C.columns.get_loc('Font_Interpreter')]

    return plot_parameters

