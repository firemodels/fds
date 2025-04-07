import numpy as np
import matplotlib.pyplot as plt
import os, sys
import pandas as pd
from matplotlib.ticker import ScalarFormatter

# Suppress RankWarning specifically
import warnings
warnings.simplefilter('ignore', np.RankWarning)

# include FDS plot styles, etc.
filedir = os.path.dirname(__file__)
firemodels = os.path.join(filedir,'..','..','..','..')
sys.path.append(filedir+os.sep+'..'+os.sep)
import fdsplotlib

# Get plot style parameters
plot_style = fdsplotlib.get_plot_style("fds")

# Define paths
base_path = os.path.join(firemodels,'out','USFS_Catchpole')
fig_path = os.path.join(firemodels,'fds','Manuals','FDS_Validation_Guide','SCRIPT_FIGURES','USFS_Catchpole')
validation_path = os.path.join(firemodels,'fds','Validation','USFS_Catchpole','FDS_Input_Files')

# Experiment parameters
tests = pd.read_csv(os.path.join(validation_path,"Test_Matrix.csv"))

for ti,test in tests.iterrows():
    R = tests['R'].iloc[ti]
    chid = tests['Test'].iloc[ti]
    fds_file = os.path.join(base_path, f"{chid}_devc.csv")
    git_file = os.path.join(base_path, f"{chid}_git.txt")
    fig_file = os.path.join(fig_path, f"{chid}.pdf")

    print("  plotting {}...".format(chid))
    
    if os.path.exists(fds_file) is False:
        print(f'Error: File {fds_file} does not exist. Skipping case.')
        continue
    
    fds_data = pd.read_csv(fds_file,header=1)
    
    # fit spread rate
    try:
        # fit slope filtering to positions greater than 2 m from ignition
        R_FDS,intercept = np.polyfit(fds_data[fds_data['x']>2]['Time'], 
                                 fds_data[fds_data['x']>2]['x'], 1)  # 1 indicates linear fit (degree 1)
    # not enough data to fit 
    except:
        R_FDS=0.
        
    if R_FDS<0:
        R_FDS=0
    
    # Exp results
    x_exp = np.array([0., 8.])
    t_exp = np.array([0., 8./R])
    version_string = fdsplotlib.get_version_string(git_file)
    fig = fdsplotlib.plot_to_fig(x_data=t_exp, y_data=x_exp,
                                 data_label='EXP',
                                 marker_style='k-')
    
    fig = fdsplotlib.plot_to_fig(x_data=fds_data['Time'].values, y_data=fds_data['x'].values,
                                 revision_label=version_string,
                                 figure_handle=fig,
                                 data_label='FDS',
                                 x_label='Time (s)',
                                 y_label='Distance (m)',
                                 marker_style='k--',
                                 y_min=0.,y_max=8.,
                                 x_min=0.,x_max=8./R,
                                 legend_location="lower right",
                                 show_legend='show',
                                 )
    
    ax = fig.axes[0]
    fdsplotlib.add_version_string(ax, chid, plot_type='linear', scale_x=0.05, scale_y=0.9, font_size=plot_style['Label_Font_Size'])
    
    #fig.supertitle(chid)
    fig.savefig(fig_file, backend='pdf')
    plt.close()
        
    # write table for dataplot
    tests.loc[tests['Test'].index[ti],'R_FDS'] = R_FDS
    test = tests.iloc[ti]
    test = test.drop('Test')
    out_file = os.path.join(base_path,f"{chid}_FDS.csv")
    pd.DataFrame([test]).to_csv(out_file,index=False)
    
    # add fds data to full table for summary plotting      
    tests.loc[ti,'R_FDS'] = R_FDS
    
##### Create summary plots

# variables of interest
dep_variables={"s":"Surface-to-Volume Ratio (1/m)",
           "beta":"Packing Ratio (-)",
           "U":"Wind Speed (m/s)",
           "M":"FMC (-)"}

# fuel labels for filtering data
fuel_labels=["MF","EXSC","PPMC","EX"]
for dvar in dep_variables:
    plt.rcParams['mathtext.fontset'] = 'stix' #'dejavuserif'
    plt.rcParams['mathtext.default'] = 'rm'
    plt.rcParams['axes.formatter.use_mathtext'] = False
    fig_file = os.path.join(fig_path, f"Catchpole_R_v_{dvar}.pdf")    
    
    # show +/- 20% relative error
    [xmin,xmax] = [tests[dvar].min(),tests[dvar].max()]
    colors = ['b','g','r','c','m','y','k'] # matlab defauls
    fig = fdsplotlib.plot_to_fig(x_data=[xmin, xmax], y_data=[0.8,0.8],
                             plot_type='semilogy', marker_style='k--',)
    for i in range(0, len(fuel_labels)):
        fuel = fuel_labels[i]
        filtered_data = tests[tests['Test'].str.startswith(fuel)]
        color = colors[i]
        if fuel=='EX':
            filtered_data = tests[
                (tests['Test'].str.startswith(fuel))&(~tests['Test'].str.startswith('EXSC'))]
        fig = fdsplotlib.plot_to_fig(x_data=filtered_data[dvar], y_data=filtered_data['R_FDS']/filtered_data['R'],
                                 data_label=fuel, plot_type='semilogy', marker_style=colors[i]+'o', figure_handle=fig)
    fig = fdsplotlib.plot_to_fig(x_data=[xmin, xmax], y_data=[1.2,1.2],
                             plot_type='semilogy', marker_style='k--', figure_handle=fig,
                             x_min=xmin, x_max=xmax, y_min=3e-3, y_max=1.1e1,
                             x_label=dep_variables[dvar],y_label=r"$\mathrm{R_{FDS}/R_{Exp}}$ $\mathrm{(-)}$",
                             show_legend='show', legend_framealpha=1.0, legend_location=2)
    plt.gca().yaxis.set_major_formatter(ScalarFormatter()) 
    
    # add version sting
    version_str = fdsplotlib.get_version_string(git_file)
    fdsplotlib.add_version_string(ax, version_str, plot_type='semilogy')        
    
    plt.tight_layout()
    plt.savefig(fig_file)
    plt.close()





# plot no-spread conditions

fig_file = os.path.join(fig_path, "Catchpole_no_spread.pdf")
fig, ax = plt.subplots(figsize=(plot_style["Paper_Width"], plot_style["Paper_Height"]))

# dummy column for labeling
tests['category']='go'
tests.loc[tests['R_FDS']<1e-5,'category'] = 'no-go'

# normalize by max and min
tests_normalized = tests
tests_normalized[list(dep_variables.keys())] = tests[list(dep_variables.keys())].apply(
    lambda x: (x - x.min()) / (x.max() - x.min()))

# move M toward the middle of x-axis for more clarity
tests_normalized=tests_normalized[['category','s','beta','M','U']]

pd.plotting.parallel_coordinates(tests_normalized, 'category', 
                                 cols=['s','beta','M','U'],
                                 color=[(1.,0.,0.,1), (0.,0.,0.,.2)],
                                 ax=ax,
                                 ls='-')

ax.set_ylim([0, 1])
ax.set_yticks([0, 1],['min','max'])
plt.legend(loc="upper left", fontsize=plot_style["Key_Font_Size"], 
           framealpha=1,frameon=True)

# Show the plot
plt.tight_layout()
plt.savefig(fig_file)
plt.close()

