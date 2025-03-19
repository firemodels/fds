import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

# include FDS plot styles, etc.
import fdsplotlib

# Get plot style parameters
plot_style = fdsplotlib.get_plot_style("fds")
plt.rcParams['text.usetex'] = True # supports latex math (set per plot below)
plt.rcParams["pdf.use14corefonts"] = True # forces matplotlib to write native pdf fonts rather than embed
plt.rcParams["font.family"] = plot_style["Font_Name"]
plt.rcParams["font.size"] = plot_style["Label_Font_Size"]

# Define paths
base_path = "../../../out/USFS_Catchpole/"
fig_path = "../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/USFS_Catchpole/"

# Experiment parameters
tests = pd.read_csv("../../Validation/USFS_Catchpole/FDS_Input_Files/Test_Matrix.csv")

for ti,test in tests.iterrows():
    chid = test['Test']
    fds_file = os.path.join(base_path, f"{chid}_devc.csv")
    git_file = os.path.join(base_path, f"{chid}_git.txt")
    fig_file = os.path.join(fig_path, f"{chid}.pdf")
    
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
    
    fig, ax = plt.subplots(figsize=(plot_style["Paper_Width"], plot_style["Paper_Height"]))
    
    # Exp results
    x_exp = np.array([0., 8.])
    t_exp = np.array([0., 8./test['R']])
    ax.plot(t_exp,x_exp,'k-',label='Exp')
    
    # FDS results
    ax.plot(fds_data['Time'],fds_data['x'],'k--',label='FDS')
    
    # plot attributes
    ax.set_xlabel("Time (s)",
                  fontdict={"fontname": plot_style["Font_Name"], 
                            "fontsize": plot_style["Label_Font_Size"]})
    ax.set_ylabel("Distance (m)",
                  fontdict={"fontname": plot_style["Font_Name"], 
                            "fontsize": plot_style["Label_Font_Size"]})
    t_end = max(fds_data['Time'].max(),8./test['R'])
    ax.set_xlim([0, t_end])
    ax.set_ylim([0, 8.])
    plt.legend(loc="lower right", fontsize=plot_style["Key_Font_Size"], 
               framealpha=1,frameon=True)
    ax.set_title(chid,fontsize=plot_style["Title_Font_Size"],
                 loc="left",x=0.05,y=0.9)

    # add version sting
    version_str = fdsplotlib.get_version_string(git_file)
    fdsplotlib.add_version_string(ax, version_str, plot_type='linear')        

    fig.tight_layout()
    plt.savefig(fig_file)
    plt.close()
        
    # write table for dataplot
    test['R_FDS'] = R_FDS
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
    
    fig_file = os.path.join(fig_path, f"Catchpole_R_v_{dvar}.pdf")    
    fig, ax = plt.subplots(figsize=(plot_style["Paper_Width"], plot_style["Paper_Height"]))
    
    # show +/- 20% relative error
    [xmin,xmax] = [tests[dvar].min(),tests[dvar].max()]
    ax.semilogy([xmin,xmax],[0.8,0.8],'k--')
    ax.semilogy([xmin,xmax],[1.2,1.2],'k--')
    
    for fuel in fuel_labels:
        filtered_data = tests[tests['Test'].str.startswith(fuel)]
        if fuel=='EX':
            filtered_data = tests[
                (tests['Test'].str.startswith(fuel))&(~tests['Test'].str.startswith('EXSC'))]
    
        ax.semilogy(filtered_data[dvar],filtered_data['R_FDS']/filtered_data['R'],
                    '.',label=fuel)
    
    # plot attributes
    ax.set_xlabel(dep_variables[dvar],
                  fontdict={"fontname": plot_style["Font_Name"], 
                            "fontsize": plot_style["Label_Font_Size"]})
    ax.set_ylabel("$R_{FDS}/R_{Exp}$ (-)",
                  fontdict={"fontname": plot_style["Font_Name"], 
                            "fontsize": plot_style["Label_Font_Size"]})
    plt.legend(fontsize=plot_style["Key_Font_Size"], 
               framealpha=1,frameon=True)
    ax.set_xlim([xmin,xmax])
    
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

