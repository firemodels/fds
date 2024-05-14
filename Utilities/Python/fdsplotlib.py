#!/usr/bin/env python
"""
fdsplotlib.py
by Randy McDermott
April 2024

Fire Dynamics Simulator (FDS) Plot Library

Collection of functions for plotting and analysis
"""

import sys
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns # colorblind pallete

def plot_to_fig(x_data,y_data,**kwargs):
    """
    Create a simple x,y plot and return the fig handle
    """
    # # useful debug statements
    # print(x_data)
    # print(y_data)
    # for key, value in kwargs.items():
    #     print ("%s == %s" %(key, value))

    ##### default parameters ######
    default_figure_size = (8,6)
    default_ticklabel_fontsize = 16
    default_axeslabel_fontsize = 18
    default_legend_fontsize = 16
    default_title_fontsize = 18
    default_subtitle_fontsize = 16
    default_stamp_fontsize = 12
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

    if kwargs.get('plot_type'):
        plot_type=kwargs.get('plot_type')
    else:
        plot_type='linear'

    # generate the main x,y plot
    if plot_type=='linear':
        ax.plot(x_data,y_data,
            markevery=kwargs.get('data_markevery'),
            label=kwargs.get('data_label'),
            markerfacecolor=kwargs.get('marker_fill_color'),
            markeredgecolor=kwargs.get('marker_edge_color'),
            markeredgewidth=kwargs.get('line_width'),
            marker=kwargs.get('marker_style'),
            markersize=kwargs.get('marker_size'),
            linestyle=kwargs.get('line_style'),
            linewidth=kwargs.get('line_width'),
            color=kwargs.get('line_color'))

    if plot_type=='loglog':
        ax.loglog(x_data,y_data,
            markevery=kwargs.get('data_markevery'),
            label=kwargs.get('data_label'),
            markerfacecolor=kwargs.get('marker_fill_color'),
            markeredgecolor=kwargs.get('marker_edge_color'),
            markeredgewidth=kwargs.get('line_width'),
            marker=kwargs.get('marker_style'),
            markersize=kwargs.get('marker_size'),
            linestyle=kwargs.get('line_style'),
            linewidth=kwargs.get('line_width'),
            color=kwargs.get('line_color'))

    if plot_type=='semilogx':
        ax.semilogx(x_data,y_data,
            markevery=kwargs.get('data_markevery'),
            label=kwargs.get('data_label'),
            markerfacecolor=kwargs.get('marker_fill_color'),
            markeredgecolor=kwargs.get('marker_edge_color'),
            markeredgewidth=kwargs.get('line_width'),
            marker=kwargs.get('marker_style'),
            markersize=kwargs.get('marker_size'),
            linestyle=kwargs.get('line_style'),
            linewidth=kwargs.get('line_width'),
            color=kwargs.get('line_color'))

    if plot_type=='semilogy':
        ax.semilogy(x_data,y_data,
            markevery=kwargs.get('data_markevery'),
            label=kwargs.get('data_label'),
            markerfacecolor=kwargs.get('marker_fill_color'),
            markeredgecolor=kwargs.get('marker_edge_color'),
            markeredgewidth=kwargs.get('line_width'),
            marker=kwargs.get('marker_style'),
            markersize=kwargs.get('marker_size'),
            linestyle=kwargs.get('line_style'),
            linewidth=kwargs.get('line_width'),
            color=kwargs.get('line_color'))

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

    if kwargs.get('subtitle_fontsize'):
        subtitle_fontsize=kwargs.get('subtitle_fontsize')
    else:
        subtitle_fontsize=default_subtitle_fontsize

    # set axes and tick properties
    if kwargs.get('x_min')==None or kwargs.get('x_max')==None or math.isnan(kwargs.get('x_min')) or math.isnan(kwargs.get('x_max')):
        xmin=min(x_data) - 0.05*(max(x_data)-min(x_data))
        xmax=max(x_data) + 0.05*(max(x_data)-min(x_data))
    else:
        xmin=kwargs.get('x_min')
        xmax=kwargs.get('x_max')

    if kwargs.get('y_min')==None or kwargs.get('y_max')==None or math.isnan(kwargs.get('y_min')) or math.isnan(kwargs.get('y_max')):
        ymin=min(y_data) - 0.05*(max(y_data)-min(y_data))
        ymax=max(y_data) + 0.05*(max(y_data)-min(y_data))
    else:
        ymin=kwargs.get('y_min')
        ymax=kwargs.get('y_max')

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    # if axis labels are too large or small, use sci notation

    if plot_type=='linear':
        x_axis_exponent_min = -3
        x_axis_exponent_max = 3
        if kwargs.get('figure_x_axis_exponent_min'):
            x_axis_exponent_min=kwargs.get('figure_x_axis_exponent_min')
        if kwargs.get('figure_x_axis_exponent_max'):
            x_axis_exponent_max=kwargs.get('figure_x_axis_exponent_max')
        ax.ticklabel_format(axis='x',scilimits=(x_axis_exponent_min,x_axis_exponent_max))

    if plot_type=='linear':
        y_axis_exponent_min = -3
        y_axis_exponent_max = 3
        if kwargs.get('figure_y_axis_exponent_min'):
            y_axis_exponent_min=kwargs.get('figure_y_axis_exponent_min')
        if kwargs.get('figure_y_axis_exponent_max'):
            y_axis_exponent_max=kwargs.get('figure_y_axis_exponent_max')
        ax.ticklabel_format(axis='y',scilimits=(y_axis_exponent_min,y_axis_exponent_max))

    if kwargs.get('x_nticks'):
        ax.set_xticks(np.linspace(start = kwargs.get('x_min'), stop = kwargs.get('x_max'), num = kwargs.get('x_nticks'), endpoint=True))

    if kwargs.get('y_nticks'):
        ax.set_yticks(np.linspace(start = kwargs.get('y_min'), stop = kwargs.get('y_max'), num = kwargs.get('y_nticks'), endpoint=True))

    # if plot_type=='linear':
    #     xpos=xmin+0.05*(xmax-xmin)
    #     ypos1=ymin+0.900*(ymax-ymin)
    #     ypos2=ymin+0.825*(ymax-ymin)
    #     ax.text(xpos,ypos1, kwargs.get('plot_title'), fontsize=title_fontsize)
    #     ax.text(xpos,ypos2, kwargs.get('plot_subtitle'), fontsize=subtitle_fontsize)

    #     # plot Institute + MaCFP stamp
    #     ax.text(xmin+0.025*(xmax-xmin),ymax+0.01*(ymax-ymin), kwargs.get('institute_label'), fontsize=default_stamp_fontsize)
    #     ax.text(xmax-0.025*(xmax-xmin),ymax+0.01*(ymax-ymin), 'MaCFP 2021', fontsize=default_stamp_fontsize, ha='right')

    # this relative method seems to work for all plot types (linear,loglog,etc.)
    xpos=0.05
    ypos1=0.900
    ypos2=0.825
    ax.text(xpos,ypos1, kwargs.get('plot_title'), fontsize=title_fontsize, transform=ax.transAxes)
    ax.text(xpos,ypos2, kwargs.get('plot_subtitle'), fontsize=subtitle_fontsize, transform=ax.transAxes)

    # plot Institute + MaCFP (or revision) stamp
    inst_label_x = 0.025
    if len(str(ymax))>abs(y_axis_exponent_max) or len(str(ymin))>abs(y_axis_exponent_min):
        inst_label_x = 0.06 # else the institute label overlays the exponential notation multiplier
    ax.text(inst_label_x,1.01, kwargs.get('institute_label'), fontsize=default_stamp_fontsize, transform=ax.transAxes)
    ax.text(0.975,       1.01, kwargs.get('revision_label'), fontsize=default_stamp_fontsize, ha='right', transform=ax.transAxes)

    # note: this absolute method works better than fig.tight_layout(), which may change for each call of the figure
    left_adjust, bottom_adjust = get_subplots_adjust_parameters(ticklabel_fontsize,axeslabel_fontsize)
    if ymin<0:
        left_adjust = left_adjust + 0.025 # account for negative sign

    if kwargs.get('figure_left_adjust'):
        left_adjust=kwargs.get('figure_left_adjust')

    if kwargs.get('figure_bottom_adjust'):
        bottom_adjust=kwargs.get('figure_bottom_adjust')

    if kwargs.get('figure_right_adjust'):
        right_adjust=kwargs.get('figure_right_adjust')
    else:
        right_adjust = 0.95

    if kwargs.get('figure_top_adjust'):
        top_adjust=kwargs.get('figure_top_adjust')
    else:
        top_adjust = 0.95

    fig.subplots_adjust(left=left_adjust, bottom=bottom_adjust, right=right_adjust, top=top_adjust, wspace=0.2, hspace=0.2)

    return fig


def get_subplots_adjust_parameters(ticklabel_fontsize,axeslabel_fontsize):
    """
    automatically select axes adjustments based on fontsize,
    this is a brute force method, but I could not find a more
    automatic way to do it since tight_layout() is not working
    """
    left_adjust = 0.125
    bottom_adjust = 0.125

    if ticklabel_fontsize + axeslabel_fontsize > 34:
        # determined this adjustment by trial and error
        left_adjust   = 0.125 + 0.006875*(ticklabel_fontsize + axeslabel_fontsize - 34)
        bottom_adjust = 0.125 + 0.001875*(ticklabel_fontsize + axeslabel_fontsize - 34)

    return left_adjust, bottom_adjust



def get_nticks(x_min,x_max,x_tick):
    """
    converts float x_tick, if it exists, to integer nticks
    """
    nticks = int((float(x_max)-float(x_min))/float(x_tick)) + 1

    return nticks


def define_plot_parameters(C,irow):
    """
    gathers parameters from config file
    """
    class plot_parameters:

        def __repr__(self):
            return str(self)

        try:
            Dataname       = C.values[irow,C.columns.get_loc('Dataname')]
        except:
            Dataname       = 'null'

        try:
            VerStr_Filename       = C.values[irow,C.columns.get_loc('VerStr_Filename')]
        except:
            VerStr_Filename       = None

        try:
            Exp_Filename          = C.values[irow,C.columns.get_loc('Exp_Filename')]
        except:
            sys.exit('Required column header missing: Exp_Filename')

        try:
            Exp_Header_Row        = int(C.values[irow,C.columns.get_loc('Exp_Header_Row')])
        except:
            Exp_Header_Row        = 1

        try:
            Exp_Data_Row          = int(C.values[irow,C.columns.get_loc('Exp_Data_Row')])
        except:
            Exp_Data_Row          = Exp_Header_Row+1

        try:
            Exp_x_Col_Name        = C.values[irow,C.columns.get_loc('Exp_x_Col_Name')]
        except:
            sys.exit('Required column header missing: Exp_x_Col_Name')

        try:
            Exp_y_Col_Name        = C.values[irow,C.columns.get_loc('Exp_y_Col_Name')]
        except:
            sys.exit('Required column header missing: Exp_y_Col_Name')

        try:
            Exp_Data_Markevery    = C.values[irow,C.columns.get_loc('Exp_Data_Markevery')]
            if pd.isnull(Exp_Data_Markevery):
                Exp_Data_Markevery = 1
        except:
            Exp_Data_Markevery    = 1

        try:
            Exp_Data_Label        = C.values[irow,C.columns.get_loc('Exp_Data_Label')]
        except:
            Exp_Data_Label        = 'Exp'

        # what's the difference between 'None' and None?
        # 'None' is a string that will not print a line in the call to axes.plot()
        # None is a type that let's axes.plot use its default
        # None gets filled in by df.fillna() for None values in existing columns in the DataFrame
        # Therefore, in that case, the user probably means they do not want a value
        # The except block is the case where the column does not exist
        # Here we assume the user wants axes.plot to use its default

        try:
            Exp_Marker_Style = C.values[irow,C.columns.get_loc('Exp_Marker_Style')]
            if Exp_Marker_Style==None or pd.isnull(Exp_Marker_Style):
                Exp_Marker_Style = 'None'
        except:
            Exp_Marker_Style = 'o'

        try:
            Exp_Marker_Edge_Color = C.values[irow,C.columns.get_loc('Exp_Marker_Edge_Color')]
            if pd.isnull(Exp_Marker_Edge_Color):
                Exp_Marker_Edge_Color = 'black'
        except:
            Exp_Marker_Edge_Color = 'black'

        try:
            Exp_Marker_Fill_Color = C.values[irow,C.columns.get_loc('Exp_Marker_Fill_Color')]
            if pd.isnull(Exp_Marker_Fill_Color):
                Exp_Marker_Fill_Color = 'None'
        except:
            Exp_Marker_Fill_Color = 'None'

        try:
            Exp_Marker_Size = C.values[irow,C.columns.get_loc('Exp_Marker_Size')]
            if pd.isnull(Exp_Marker_Size):
                Exp_Marker_Size = 6
        except:
            Exp_Marker_Size = 6

        try:
            Exp_Line_Style = C.values[irow,C.columns.get_loc('Exp_Line_Style')]
            if Exp_Line_Style==None or pd.isnull(Exp_Line_Style):
                Exp_Line_Style = 'None'
        except:
            # if no column exists, we assume just markers for experimental data
            Exp_Line_Style = 'None'

        try:
            Exp_Line_Color = C.values[irow,C.columns.get_loc('Exp_Line_Color')]
            if Exp_Line_Color==None or pd.isnull(Exp_Line_Color):
                Exp_Line_Color = 'None'
        except:
            if Exp_Line_Style:
                Exp_Line_Color = 'black'
            else:
                Exp_Line_Color = 'None'

        try:
            Exp_Line_Width = C.values[irow,C.columns.get_loc('Exp_Line_Width')]
            if Exp_Line_Width==None or pd.isnull(Exp_Line_Width):
                Exp_Line_Width = 1.
        except:
            Exp_Line_Width = 1.

        try:
            Exp_Error_Absolute = C.values[irow,C.columns.get_loc('Exp_Error_Absolute')]
            if pd.isnull(Exp_Error_Absolute):
                Exp_Error_Absolute = 0.
        except:
            Exp_Error_Absolute = 0.

        try:
            Exp_Error_Relative = C.values[irow,C.columns.get_loc('Exp_Error_Relative')]
            if pd.isnull(Exp_Error_Relative):
                Exp_Error_Relative = 0.
        except:
            Exp_Error_Relative = 0.

        try:
            Exp_y_Error_Col_Name = C.values[irow,C.columns.get_loc('Exp_y_Error_Col_Name')]
            if math.isnan(Exp_y_Error_Col_Name): Exp_y_Error_Col_Name = None
        except:
            Exp_y_Error_Col_Name = None

        try:
            Cmp_Filename = C.values[irow,C.columns.get_loc('Cmp_Filename')]
        except:
            sys.exit('Required column header missing: Cmp_Filename')

        try:
            Cmp_Header_Row = C.values[irow,C.columns.get_loc('Cmp_Header_Row')]
            if pd.isnull(Cmp_Header_Row):
                Cmp_Header_Row = 1
        except:
            Cmp_Header_Row = 1

        try:
            Cmp_Data_Row = C.values[irow,C.columns.get_loc('Cmp_Data_Row')]
            if pd.isnull(Cmp_Data_Row):
                Cmp_Data_Row = Cmp_Header_Row+1
        except:
            Cmp_Data_Row = Cmp_Header_Row+1

        try:
            Cmp_x_Col_Name = C.values[irow,C.columns.get_loc('Cmp_x_Col_Name')]
        except:
            sys.exit('Required column header missing: Cmp_x_Col_Name')

        try:
            Cmp_y_Col_Name = C.values[irow,C.columns.get_loc('Cmp_y_Col_Name')]
        except:
            sys.exit('Required column header missing: Cmp_y_Col_Name')

        try:
            Cmp_Data_Markevery = C.values[irow,C.columns.get_loc('Cmp_Data_Markevery')]
            if pd.isnull(Cmp_Data_Markevery):
                Cmp_Data_Markevery = 1
        except:
            Cmp_Data_Markevery = 1

        try:
            Cmp_Data_Label = C.values[irow,C.columns.get_loc('Cmp_Data_Label')]
        except:
            Cmp_Data_Label = 'Cmp'

        try:
            Cmp_Marker_Style = C.values[irow,C.columns.get_loc('Cmp_Marker_Style')]
            if pd.isnull(Cmp_Marker_Style):
                Cmp_Marker_Style = 'None'
        except:
            Cmp_Marker_Style = 'None'

        try:
            Cmp_Marker_Edge_Color = C.values[irow,C.columns.get_loc('Cmp_Marker_Edge_Color')]
            if pd.isnull(Cmp_Marker_Edge_Color):
                Cmp_Marker_Edge_Color = 'black'
        except:
            Cmp_Marker_Edge_Color = 'black'

        try:
            Cmp_Marker_Fill_Color = C.values[irow,C.columns.get_loc('Cmp_Marker_Fill_Color')]
            if pd.isnull(Cmp_Marker_Fill_Color):
                Cmp_Marker_Fill_Color = 'None'
        except:
            Cmp_Marker_Fill_Color = 'None'

        try:
            Cmp_Marker_Size = C.values[irow,C.columns.get_loc('Cmp_Marker_Size')]
            if pd.isnull(Cmp_Marker_Size):
                Cmp_Marker_Size = 6
        except:
            Cmp_Marker_Size = 6

        try:
            Cmp_Line_Style = C.values[irow,C.columns.get_loc('Cmp_Line_Style')]
            if Cmp_Line_Style==None or pd.isnull(Cmp_Line_Style):
                Cmp_Line_Style = 'None'
        except:
            Cmp_Line_Style = 'None'

        try:
            Cmp_Line_Color = C.values[irow,C.columns.get_loc('Cmp_Line_Color')]
            if pd.isnull(Cmp_Line_Color):
                Cmp_Line_Color = 'black'
        except:
            Cmp_Line_Color = 'black'

        try:
            Cmp_Line_Width = C.values[irow,C.columns.get_loc('Cmp_Line_Width')]
            if Cmp_Line_Width==None or pd.isnull(Cmp_Line_Width):
                Cmp_Line_Width = 1.
        except:
            Cmp_Line_Width = 1.

        try:
            Plot_Title = C.values[irow,C.columns.get_loc('Plot_Title')]
            if pd.isnull(Plot_Title):
                Plot_Title = ''
        except:
            Plot_Title = ''

        try:
            Plot_Subtitle = C.values[irow,C.columns.get_loc('Plot_Subtitle')]
            if pd.isnull(Plot_Subtitle):
                Plot_Subtitle = ''
        except:
            Plot_Subtitle = ''

        try:
            Plot_x_Label = C.values[irow,C.columns.get_loc('Plot_x_Label')]
            if pd.isnull(Plot_x_Label):
                Plot_x_Label = 'x'
        except:
            Plot_x_Label = 'x'
        try:
            Plot_y_Label = C.values[irow,C.columns.get_loc('Plot_y_Label')]
            if pd.isnull(Plot_y_Label):
                Plot_y_Label = 'y'
        except:
            Plot_y_Label = 'y'

        try:
            Plot_x_Min = C.values[irow,C.columns.get_loc('Plot_x_Min')]
            if pd.isnull(Plot_x_Min):
                Plot_x_Min = None
        except:
            Plot_x_Min = None

        try:
            Plot_x_Max = C.values[irow,C.columns.get_loc('Plot_x_Max')]
            if pd.isnull(Plot_x_Max):
                Plot_x_Max = None
        except:
            Plot_x_Max = None

        try:
            Plot_x_Tick = C.values[irow,C.columns.get_loc('Plot_x_Tick')]
            if pd.isnull(Plot_x_Tick):
                Plot_x_Tick = None
        except:
            Plot_x_Tick = None

        try:
            Plot_y_Min = C.values[irow,C.columns.get_loc('Plot_y_Min')]
            if pd.isnull(Plot_y_Min):
                Plot_y_Min = None
        except:
            Plot_y_Min = None

        try:
            Plot_y_Max = C.values[irow,C.columns.get_loc('Plot_y_Max')]
            if pd.isnull(Plot_y_Max):
                Plot_y_Max = None
        except:
            Plot_y_Max = None

        try:
            Plot_y_Tick = C.values[irow,C.columns.get_loc('Plot_y_Tick')]
            if pd.isnull(Plot_y_Tick):
                Plot_y_Tick = None
        except:
            Plot_y_Tick = None

        try:
            Plot_x_Nticks = C.values[irow,C.columns.get_loc('Plot_x_Nticks')]
            if pd.isnull(Plot_x_Nticks):
                Plot_x_Nticks = None
        except:
            Plot_x_Nticks = None

        try:
            if float(Plot_x_Tick)>0. and float(Plot_x_Tick)<1.e10 and float(Plot_x_Min)>-1.e10 and float(Plot_x_Max)<1.e10:
                Plot_x_Nticks = get_nticks(Plot_x_Min,Plot_x_Max,Plot_x_Tick)
            if pd.isnull(Plot_x_Nticks):
                Plot_x_Nticks = None
        except:
            Plot_x_Nticks = None

        try:
            Plot_y_Nticks = C.values[irow,C.columns.get_loc('Plot_y_Nticks')]
            if pd.isnull(Plot_y_Nticks):
                Plot_y_Nticks = None
        except:
            Plot_y_Nticks = None

        try:
            if float(Plot_y_Tick)>0. and float(Plot_y_Tick)<1.e10 and float(Plot_y_Min)>-1.e10 and float(Plot_y_Max)<1.e10:
                Plot_y_Nticks = get_nticks(Plot_y_Min,Plot_y_Max,Plot_y_Tick)
            if pd.isnull(Plot_y_Nticks):
                Plot_y_Nticks = None
        except:
            Plot_y_Nticks = None

        try:
            Plot_x_Axis_Exponent_Min = C.values[irow,C.columns.get_loc('Plot_x_Axis_Exponent_Min')]
            if pd.isnull(Plot_x_Axis_Exponent_Min):
                Plot_x_Axis_Exponent_Min = None
        except:
            Plot_x_Axis_Exponent_Min = None

        try:
            Plot_x_Axis_Exponent_Max = C.values[irow,C.columns.get_loc('Plot_x_Axis_Exponent_Max')]
            if pd.isnull(Plot_x_Axis_Exponent_Max):
                Plot_x_Axis_Exponent_Max = None
        except:
            Plot_x_Axis_Exponent_Max = None

        try:
            Plot_y_Axis_Exponent_Min = C.values[irow,C.columns.get_loc('Plot_y_Axis_Exponent_Min')]
            if pd.isnull(Plot_y_Axis_Exponent_Min):
                Plot_y_Axis_Exponent_Min = None
        except:
            Plot_y_Axis_Exponent_Min = None

        try:
            Plot_y_Axis_Exponent_Max = C.values[irow,C.columns.get_loc('Plot_y_Axis_Exponent_Max')]
            if pd.isnull(Plot_y_Axis_Exponent_Max):
                Plot_y_Axis_Exponent_Max = None
        except:
            Plot_y_Axis_Exponent_Max = None

        try:
            Plot_Flip_Axis = C.values[irow,C.columns.get_loc('Plot_Flip_Axis')]
            if pd.isnull(Plot_Flip_Axis):
                Plot_Flip_Axis = False
        except:
            Plot_Flip_Axis = False

        try:
            Plot_Show_Legend = C.values[irow,C.columns.get_loc('Plot_Show_Legend')]
            if pd.isnull(Plot_Show_Legend):
                Plot_Show_Legend = True
        except:
            Plot_Show_Legend = True

        try:
            Plot_Legend_Location = C.values[irow,C.columns.get_loc('Plot_Legend_Location')]
            if Plot_Legend_Location==None or pd.isnull(Plot_Legend_Location):
                Plot_Legend_Location = 'best'
        except:
            Plot_Legend_Location = 'best'

        if Plot_Legend_Location.isdigit():
            Plot_Legend_Location=int(Plot_Legend_Location)

        try:
            Plot_Filename = C.values[irow,C.columns.get_loc('Plot_Filename')]
        except:
            sys.exit('Required column header missing: Plot_Filename')

        try:
            Plot_Left_Adjust = C.values[irow,C.columns.get_loc('Plot_Left_Adjust')]
            if pd.isnull(Plot_Left_Adjust):
                Plot_Left_Adjust = None
        except:
            Plot_Left_Adjust = None

        try:
            Plot_Right_Adjust = C.values[irow,C.columns.get_loc('Plot_Right_Adjust')]
            if pd.isnull(Plot_Right_Adjust):
                Plot_Right_Adjust = None
        except:
            Plot_Right_Adjust = None

        try:
            Plot_Bottom_Adjust = C.values[irow,C.columns.get_loc('Plot_Bottom_Adjust')]
            if pd.isnull(Plot_Bottom_Adjust):
                Plot_Bottom_Adjust = None
        except:
            Plot_Bottom_Adjust = None

        try:
            Plot_Top_Adjust = C.values[irow,C.columns.get_loc('Plot_Top_Adjust')]
            if pd.isnull(Plot_Top_Adjust):
                Plot_Top_Adjust = None
        except:
            Plot_Top_Adjust = None

        try:
            Plot_Figure_Width = C.values[irow,C.columns.get_loc('Plot_Figure_Width')]
            if Plot_Figure_Width==None or pd.isnull(Plot_Figure_Width):
                Plot_Figure_Width=8
        except:
            Plot_Figure_Width = 8

        try:
            Plot_Figure_Height = C.values[irow,C.columns.get_loc('Plot_Figure_Height')]
            if Plot_Figure_Height==None or pd.isnull(Plot_Figure_Height):
                Plot_Figure_Height=6
        except:
            Plot_Figure_Height = 6

        # try:
        #     Seaborn_Color = C.values[irow,C.columns.get_loc('Seaborn_Color')]
        #     if pd.isnull(Seaborn_Color):
        #         Seaborn_Color = False
        # except:
        #     Seaborn_Color = False

    return plot_parameters


def dataplot(config_filename,**kwargs):

    import os
    import matplotlib.pyplot as plt
    import pandas as pd
    # from matplotlib import rc
    # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    # rc('font',**{'family':'serif','serif':['Palatino']})

    plt.rcParams['mathtext.fontset'] = 'dejavusans' # valid strings are ['dejavusans', 'dejavuserif', 'cm', 'stix', 'stixsans', 'custom']
    plt.rcParams['pdf.compression'] = 6 # 0 to 9 (file size was no smaller with 9)
    plt.rcParams['pdf.fonttype'] = 3 # 42 embeds true type fonts into pdf (large file size)
    # plt.rcParams['text.usetex'] = True # supports latex math
    plt.rcParams["pdf.use14corefonts"] = True # forces matplotlib to write native pdf fonts rather then embed

    # defaults
    institute = ''
    revision  = ''
    configdir = ''
    expdir = ''
    cmpdir = ''
    pltdir = ''
    close_figs = False
    verbose = False
    plot_list = ['all']
    plot_range = range(10000)

    if kwargs.get('institute'):
        institute = kwargs.get('institute')

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
        plot_range = range(plot_range_in[0]-2,plot_range_in[-1])

    # read the config file
    df = pd.read_csv(configdir+config_filename, sep=' *, *', engine='python', comment='#', quotechar='"')
    C = df.where(pd.notnull(df), None)

    Plot_Filename_Last = None
    Exp_Data_Label_Last = None
    f_Last = plt.figure()

    # loop over the rows of the config file
    for irow in C.index:

        if irow not in plot_range:
            continue

        try:

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
                E = pd.read_csv(expdir+pp.Exp_Filename, header=pp.Exp_Header_Row-1, sep=' *, *', engine='python')
                if (pp.Exp_Data_Row-pp.Exp_Header_Row==1):
                    x = E[pp.Exp_x_Col_Name].values[:].astype(float)
                    y = E[pp.Exp_y_Col_Name].values[:].astype(float)
                    if (pp.Exp_y_Error_Col_Name):
                        y_error = E[pp.Exp_y_Error_Col_Name].values[:].astype(float)
                    else:
                        y_error = 0.
                else:
                    # don't exactly understand this, but df.values behave differently if they are object type
                    # when the header and data rows are separated, then there are usually strings in the df.values
                    x = E[pp.Exp_x_Col_Name].values[pp.Exp_Data_Row-2:].astype(float)
                    y = E[pp.Exp_y_Col_Name].values[pp.Exp_Data_Row-2:].astype(float)
                    if (pp.Exp_y_Error_Col_Name):
                        y_error = E[pp.Exp_y_Error_Col_Name].values[pp.Exp_Data_Row-2:].astype(float)
                    else:
                        y_error = 0.

                # plot the exp data
                f = plot_to_fig(x_data=x, y_data=y,
                    data_markevery=pp.Exp_Data_Markevery,
                    data_label=pp.Exp_Data_Label,
                    y_error_absolute=pp.Exp_Error_Absolute,
                    y_error_relative=pp.Exp_Error_Relative,
                    y_error_vector=y_error,
                    x_label=pp.Plot_x_Label,
                    y_label=pp.Plot_y_Label,
                    marker_style=pp.Exp_Marker_Style,
                    marker_fill_color=pp.Exp_Marker_Fill_Color,
                    marker_edge_color=pp.Exp_Marker_Edge_Color,
                    marker_size=pp.Exp_Marker_Size,
                    line_style=pp.Exp_Line_Style,
                    line_color=pp.Exp_Line_Color,
                    line_width=pp.Exp_Line_Width,
                    x_min=pp.Plot_x_Min,x_max=pp.Plot_x_Max,x_nticks=pp.Plot_x_Nticks,
                    y_min=pp.Plot_y_Min,y_max=pp.Plot_y_Max,y_nticks=pp.Plot_y_Nticks,
                    show_legend=pp.Plot_Show_Legend,legend_location=pp.Plot_Legend_Location,
                    figure_size=(pp.Plot_Figure_Width,pp.Plot_Figure_Height),
                    figure_left_adjust=pp.Plot_Left_Adjust,
                    figure_right_adjust=pp.Plot_Right_Adjust,
                    figure_bottom_adjust=pp.Plot_Bottom_Adjust,
                    figure_top_adjust=pp.Plot_Top_Adjust,
                    figure_x_axis_exponent_min=pp.Plot_x_Axis_Exponent_Min,
                    figure_x_axis_exponent_max=pp.Plot_x_Axis_Exponent_Max,
                    figure_y_axis_exponent_min=pp.Plot_y_Axis_Exponent_Min,
                    figure_y_axis_exponent_max=pp.Plot_y_Axis_Exponent_Max
                    )

                # plt.figure(f.number) # make figure current
                # plt.show()
            else:
                f = f_Last

                if pp.Exp_Data_Label!=Exp_Data_Label_Last:

                    # read data from exp file
                    # set header to the row where column names are stored (Python is 0 based)
                    E = pd.read_csv(expdir+pp.Exp_Filename, header=pp.Exp_Header_Row-1, sep=' *, *', engine='python')

                    if (pp.Exp_Data_Row-pp.Exp_Header_Row==1):
                        x = E[pp.Exp_x_Col_Name].values[:].astype(float)
                        y = E[pp.Exp_y_Col_Name].values[:].astype(float)
                        if (pp.Exp_y_Error_Col_Name):
                            y_error = E[pp.Exp_y_Error_Col_Name].values[:].astype(float)
                        else:
                            y_error = 0.
                    else:
                        # don't exactly understand this, but df.values behave differently if they are object type
                        # when the header and data rows are separated, then there are usually strings in the df.values
                        x = E[pp.Exp_x_Col_Name].values[pp.Exp_Data_Row-2:].astype(float)
                        y = E[pp.Exp_y_Col_Name].values[pp.Exp_Data_Row-2:].astype(float)
                        if (pp.Exp_y_Error_Col_Name):
                            y_error = E[pp.Exp_y_Error_Col_Name].values[pp.Exp_Data_Row-2:].astype(float)
                        else:
                            y_error = 0.

                    # plot the exp data
                    f = plot_to_fig(x_data=x, y_data=y,
                        figure_handle=f,
                        data_markevery=pp.Exp_Data_Markevery,
                        data_label=pp.Exp_Data_Label,
                        y_error_absolute=pp.Exp_Error_Absolute,
                        y_error_relative=pp.Exp_Error_Relative,
                        y_error_vector=y_error,
                        x_label=pp.Plot_x_Label,
                        y_label=pp.Plot_y_Label,
                        marker_style=pp.Exp_Marker_Style,
                        marker_fill_color=pp.Exp_Marker_Fill_Color,
                        marker_edge_color=pp.Exp_Marker_Edge_Color,
                        marker_size=pp.Exp_Marker_Size,
                        line_style=pp.Exp_Line_Style,
                        line_color=pp.Exp_Line_Color,
                        line_width=pp.Exp_Line_Width,
                        x_min=pp.Plot_x_Min,x_max=pp.Plot_x_Max,x_nticks=pp.Plot_x_Nticks,
                        y_min=pp.Plot_y_Min,y_max=pp.Plot_y_Max,y_nticks=pp.Plot_y_Nticks,
                        show_legend=pp.Plot_Show_Legend,legend_location=pp.Plot_Legend_Location,
                        figure_size=(pp.Plot_Figure_Width,pp.Plot_Figure_Height),
                        figure_left_adjust=pp.Plot_Left_Adjust,
                        figure_right_adjust=pp.Plot_Right_Adjust,
                        figure_bottom_adjust=pp.Plot_Bottom_Adjust,
                        figure_top_adjust=pp.Plot_Top_Adjust,
                        figure_x_axis_exponent_min=pp.Plot_x_Axis_Exponent_Min,
                        figure_x_axis_exponent_max=pp.Plot_x_Axis_Exponent_Max,
                        figure_y_axis_exponent_min=pp.Plot_y_Axis_Exponent_Min,
                        figure_y_axis_exponent_max=pp.Plot_y_Axis_Exponent_Max
                        )

            # get the model results
            M = pd.read_csv(cmpdir+pp.Cmp_Filename, header=pp.Cmp_Header_Row-1, sep=' *, *', engine='python')

            if (pp.Cmp_Data_Row-pp.Cmp_Header_Row==1):
                x = M[pp.Cmp_x_Col_Name].values[:].astype(float)
                y = M[pp.Cmp_y_Col_Name].values[:].astype(float)
            else:
                x = M[pp.Cmp_x_Col_Name].values[pp.Cmp_Data_Row-2:].astype(float)
                y = M[pp.Cmp_y_Col_Name].values[pp.Cmp_Data_Row-2:].astype(float)

            version_string = revision
            if (pp.VerStr_Filename):
                file1 = open(cmpdir+pp.VerStr_Filename,"r")
                Lines = file1.readlines()
                version_string = Lines[0].strip()
                file1.close()

            f = plot_to_fig(x_data=x, y_data=y,
                institute_label=institute,
                revision_label=version_string,
                figure_handle=f,
                data_markevery=pp.Cmp_Data_Markevery,
                x_label=pp.Plot_x_Label,
                y_label=pp.Plot_y_Label,
                data_label=pp.Cmp_Data_Label,
                marker_style=pp.Cmp_Marker_Style,
                marker_fill_color=pp.Cmp_Marker_Fill_Color,
                marker_edge_color=pp.Cmp_Marker_Edge_Color,
                marker_size=pp.Cmp_Marker_Size,
                line_style=pp.Cmp_Line_Style,
                line_color=pp.Cmp_Line_Color,
                line_width=pp.Cmp_Line_Width,
                x_min=pp.Plot_x_Min,x_max=pp.Plot_x_Max,x_nticks=pp.Plot_x_Nticks,
                y_min=pp.Plot_y_Min,y_max=pp.Plot_y_Max,y_nticks=pp.Plot_y_Nticks,
                show_legend=pp.Plot_Show_Legend,legend_location=pp.Plot_Legend_Location,
                plot_title=pp.Plot_Title,
                plot_subtitle=pp.Plot_Subtitle,
                figure_left_adjust=pp.Plot_Left_Adjust,
                figure_right_adjust=pp.Plot_Right_Adjust,
                figure_bottom_adjust=pp.Plot_Bottom_Adjust,
                figure_top_adjust=pp.Plot_Top_Adjust,
                figure_x_axis_exponent_min=pp.Plot_x_Axis_Exponent_Min,
                figure_x_axis_exponent_max=pp.Plot_x_Axis_Exponent_Max,
                figure_y_axis_exponent_min=pp.Plot_y_Axis_Exponent_Min,
                figure_y_axis_exponent_max=pp.Plot_y_Axis_Exponent_Max
                )

            plt.figure(f.number) # make figure current
            # plt.show()

            # create plot directory if it does not exist
            isDir = os.path.isdir(pltdir)
            if not isDir:
                os.mkdir(pltdir)

            plt.savefig(pltdir + pp.Plot_Filename, backend='pdf')

            Plot_Filename_Last = pp.Plot_Filename
            Exp_Data_Label_Last = pp.Exp_Data_Label
            f_Last = f

        except:
            print("Error in row {whichrow}, skipping case...".format(whichrow=irow+1))
            continue

