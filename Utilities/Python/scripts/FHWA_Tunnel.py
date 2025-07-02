
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import os
import warnings
# include FDS plot styles, etc.
import fdsplotlib

warnings.filterwarnings('ignore')

# McGrattan
# 9-27-2022
# FHWA_Tunnel.py
#
# This script creates several different kinds of contour and scatter plots for the FHWA Tunnel simulations.

# clear all - not needed in Python
# close all
plt.close('all')

outdir = '../../../out/FHWA_Tunnel/'
expdir = '../../../exp/FHWA_Tunnel/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FHWA_Tunnel/'

pos = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.5, 9.5, 10.5, 11.5]
test = ['IFAB-07', 'IFAB-08', 'IFAB-09', 'IFAB-10', 'IFAB-11', 'IFAB-13', 'IFAB-14', 'IFAB-15', 'IFAB-19', 'IFAB-22', 'IFAB-24']
test2 = ['Test 7', 'Test 8', 'Test 9', 'Test 10', 'Test 11', 'Test 13', 'Test 14', 'Test 15', 'Test 19', 'Test 22', 'Test 24']
single_level = [50]
setpoint = [10000, 400, 399, 338, 240, 322, 390, 420, 360, 10000, 10000]

plot_style = fdsplotlib.get_plot_style("fds")

def addverstr(ax, git_filename, style, x, y, fontname, interpreter, fontsize):
    """Add version string to plot - placeholder implementation"""
    # This function would read git information and add it to the plot
    # For now, we'll add a placeholder text
    if os.path.exists(git_filename):
        ax.text(x, y, 'Git Version Info', transform=ax.transAxes, 
                fontname=fontname, fontsize=fontsize)

for k in range(11):  # Experiments (0-based indexing in Python)

    n_res = 1
#   if k==7:  # k==8 in MATLAB (1-based) corresponds to k==7 in Python (0-based)
#       n_res = 2

    for jj in range(n_res):

        # clear M E - handled by reassignment

        if jj == 0:
            # M = importdata([outdir,test{k},'_cat_devc.csv'],',',2);
            M_data = pd.read_csv(outdir + test[k] + '_cat_devc.csv', skiprows=2)
            M = {'data': M_data.values}
            
            # E = importdata([expdir,test{k},'_avg.csv'],',',2);
            E_data = pd.read_csv(expdir + test[k] + '_avg.csv', skiprows=2)
            E = {'data': E_data.values}
        elif jj == 1:
            # M = importdata([outdir,test{k},'_fine_cat_devc.csv'],',',2);
            M_data = pd.read_csv(outdir + test[k] + '_fine_cat_devc.csv', skiprows=2)
            M = {'data': M_data.values}
            
            # E = importdata([expdir,test{k},'_avg.csv'],',',2);
            E_data = pd.read_csv(expdir + test[k] + '_avg.csv', skiprows=2)
            E = {'data': E_data.values}

        # For each experiment, make a contour plot of the extent of a single temperature contour at each time during the experiment

        # clear X_mod Y_mod Z_mod X_exp Y_exp Z_exp

        # [X_mod,Y_mod] = meshgrid(pos(1:16),M.data(:,1)/60);
        X_mod, Y_mod = np.meshgrid(pos[0:16], M['data'][:, 0]/60)
        
        # [X_exp,Y_exp] = meshgrid(pos(1:16),E.data(:,1)/60);
        X_exp, Y_exp = np.meshgrid(pos[0:16], E['data'][:, 0]/60)
        
        # Initialize Z_mod and Z_exp
        Z_mod = np.zeros((len(M['data'][:, 0]), 16))
        Z_exp = np.zeros((len(E['data'][:, 0]), 16))
        
        for kk in range(len(M['data'][:, 0])):
            for ii in range(16):
                Z_mod[kk, ii] = M['data'][kk, ii+1]

        for kk in range(len(E['data'][:, 0])):
            for ii in range(16):
                Z_exp[kk, ii] = E['data'][kk, ii+4]

        newpoints = 100
        # [X_mod_interp,Y_mod_interp] = meshgrid(...
        #       linspace(min(min(X_mod,[],2)),max(max(X_mod,[],2)),newpoints ),...
        #       linspace(min(min(Y_mod,[],1)),max(max(Y_mod,[],1)),newpoints )...
        #     );
        X_mod_interp, Y_mod_interp = np.meshgrid(
            np.linspace(np.min(X_mod), np.max(X_mod), newpoints),
            np.linspace(np.min(Y_mod), np.max(Y_mod), newpoints)
        )
        
        # Z_mod_interp = interp2(X_mod,Y_mod,Z_mod,X_mod_interp,Y_mod_interp,'makima');
        # Flatten the arrays for griddata
        points_mod = np.column_stack((X_mod.ravel(), Y_mod.ravel()))
        values_mod = Z_mod.ravel()
        Z_mod_interp = griddata(points_mod, values_mod, 
                               (X_mod_interp, Y_mod_interp), method='cubic')

        # [X_exp_interp,Y_exp_interp] = meshgrid(...
        #       linspace(min(min(X_exp,[],2)),max(max(X_exp,[],2)),newpoints ),...
        #       linspace(min(min(Y_exp,[],1)),max(max(Y_exp,[],1)),newpoints )...
        #     );
        X_exp_interp, Y_exp_interp = np.meshgrid(
            np.linspace(np.min(X_exp), np.max(X_exp), newpoints),
            np.linspace(np.min(Y_exp), np.max(Y_exp), newpoints)
        )
        
        # Z_exp_interp = interp2(X_exp,Y_exp,Z_exp,X_exp_interp,Y_exp_interp,'makima');
        points_exp = np.column_stack((X_exp.ravel(), Y_exp.ravel()))
        values_exp = Z_exp.ravel()
        Z_exp_interp = griddata(points_exp, values_exp, 
                               (X_exp_interp, Y_exp_interp), method='cubic')

        if jj == 0:
            # reset(gca)
            # reset(gcf)
            plt.figure(figsize=(plot_style["Paper_Width"], plot_style["Paper_Height"]))
            mod_symbol = 'r-'
        else:
            mod_symbol = 'r--'

        # [C_mod,h_mod] = contour(X_mod_interp,Y_mod_interp,Z_mod_interp,single_level,mod_symbol) ; hold on
        CS_mod = plt.contour(X_mod_interp, Y_mod_interp, Z_mod_interp, single_level, colors='red', 
                            linestyles='-' if jj == 0 else '--')
        plt.clabel(CS_mod, fontsize=3, colors='red', inline_spacing=300)
        
        # [C_exp,h_exp] = contour(X_exp_interp,Y_exp_interp,Z_exp_interp,single_level,'k-') ; hold on
        CS_exp = plt.contour(X_exp_interp, Y_exp_interp, Z_exp_interp, single_level, colors='black', linestyles='-')
        plt.clabel(CS_exp, fontsize=3, colors='black', inline_spacing=300)

    # end grid resolution cases

    # plot([5.5 5.5],[0 15],'k--')
    plt.plot([5.5, 5.5], [0, 15], 'k--')
    
    # plot([0.0 15.],[setpoint(k)/60 setpoint(k)/60],'k:')
    plt.plot([0.0, 15.0], [setpoint[k]/60, setpoint[k]/60], 'k:')

    # xticks(pos)
    # xticklabels({'1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0','5.5','6.0','6.5','7.5','9.5','10.5','11.5'})
    # ax = gca;
    ax = plt.gca()
    
    # ax.XAxis.FontSize = 16;
    # ax.YAxis.FontSize = 16;
    ax.tick_params(axis='both', which='major', labelsize=16)
    
    # xlabel('Position (m)','FontSize',16,'Interpreter',Font_Interpreter)
    plt.xlabel('Position (m)', fontsize=16)
    
    # ylabel('Time (min)','FontSize',16,'Interpreter',Font_Interpreter)
    plt.ylabel('Time (min)', fontsize=16)
    
    plt.subplots_adjust(left=plot_style["Plot_X"]/plot_style["Paper_Width"], bottom=plot_style["Plot_Y"]/plot_style["Paper_Height"], 
                       right=(plot_style["Plot_X"]+plot_style["Plot_Width"])/plot_style["Paper_Width"], 
                       top=(plot_style["Plot_Y"]+plot_style["Plot_Height"])/plot_style["Paper_Height"])
    
    plt.rcParams['font.family'] = plot_style["Font_Name"]
    
    plt.axis([0, 15, 0, 15])
    
    plt.text(0.5, 4, test2[k], fontname=plot_style["Font_Name"], fontsize=10)
    plt.text(0.5, 3, 'FDS red; Exp black', fontname=plot_style["Font_Name"], fontsize=10)

    # Git_Filename = [outdir,test{k},'_cat_git.txt'];
    Git_Filename = outdir + test[k] + '_cat_git.txt'
    
    # addverstr(gca,Git_Filename,'linear',0.6,1.05,'Times','TeX',10)
    addverstr(ax, Git_Filename, 'linear', 0.6, 1.05, 'Times', 'TeX', 10)

    # set(gcf,'Units',Paper_Units);
    # set(gcf,'PaperUnits',Paper_Units);
    # set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    # set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    fig = plt.gcf()
    fig.set_size_inches(plot_style["Paper_Width"], plot_style["Paper_Height"])
    
    # print(gcf,'-dpdf',[pltdir,test{k},'_tvT'])
    plt.savefig(pltdir + test[k] + '_tvT.pdf', format='pdf', bbox_inches='tight')

    # hold off
    plt.clf()  # Clear the figure for next iteration

# end Experiment loop

print('FHWA_Tunnel completed successfully')


