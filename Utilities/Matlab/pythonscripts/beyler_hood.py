import numpy as np
import smop_util
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pgf import FigureCanvasPgf
import matplotlib as mpl
plt.rcParams['pgf.rcfonts'] = False
mpl.backend_bases.register_backend('pdf', FigureCanvasPgf)
# Floyd
# 8-28-2012
# Beyler_hood.m

plt.close('all')

expdir='../../../exp/Beyler_Hood/'
# temp/beyler_hood.m:8
outdir='../../../out/Beyler_Hood/FDS_Output_Files/'
# temp/beyler_hood.m:9
# load experimental data and FDS prediction
exp_data=np.nan_to_num(np.genfromtxt(expdir + 'Beyler_Hood_data_lr.csv', delimiter = ',')[2:,:], copy = False)
# temp/beyler_hood.m:12
Fuel=np.asarray([], dtype='object')
# temp/beyler_hood.m:14
Fuel = smop_util.safe_set(Fuel,(0,),'acetone')
# temp/beyler_hood.m:15
TestID=np.asarray([], dtype='object')
# temp/beyler_hood.m:17
TestID = smop_util.safe_set(TestID,(0,0,),117)
# temp/beyler_hood.m:18
TestID = smop_util.safe_set(TestID,(0,1,),119)
# temp/beyler_hood.m:19
TestID = smop_util.safe_set(TestID,(0,2,),122)
# temp/beyler_hood.m:20
TestID = smop_util.safe_set(TestID,(0,3,),142)
# temp/beyler_hood.m:21
TestID = smop_util.safe_set(TestID,(0,4,),145)
# temp/beyler_hood.m:22
Fuel = smop_util.safe_set(Fuel,(1,),'ethanol')
# temp/beyler_hood.m:24
TestID = smop_util.safe_set(TestID,(1,0,),106)
# temp/beyler_hood.m:26
TestID = smop_util.safe_set(TestID,(1,1,),107)
# temp/beyler_hood.m:27
TestID = smop_util.safe_set(TestID,(1,2,),108)
# temp/beyler_hood.m:28
TestID = smop_util.safe_set(TestID,(1,3,),110)
# temp/beyler_hood.m:29
TestID = smop_util.safe_set(TestID,(1,4,),115)
# temp/beyler_hood.m:30
Fuel = smop_util.safe_set(Fuel,(2,),'isopropanol')
# temp/beyler_hood.m:32
TestID = smop_util.safe_set(TestID,(2,0,),130)
# temp/beyler_hood.m:34
TestID = smop_util.safe_set(TestID,(2,1,),132)
# temp/beyler_hood.m:35
TestID = smop_util.safe_set(TestID,(2,2,),133)
# temp/beyler_hood.m:36
TestID = smop_util.safe_set(TestID,(2,3,),136)
# temp/beyler_hood.m:37
TestID = smop_util.safe_set(TestID,(2,4,),141)
# temp/beyler_hood.m:38
Fuel = smop_util.safe_set(Fuel,(3,),'methanol')
# temp/beyler_hood.m:40
TestID = smop_util.safe_set(TestID,(3,0,),942)
# temp/beyler_hood.m:42
TestID = smop_util.safe_set(TestID,(3,1,),943)
# temp/beyler_hood.m:43
TestID = smop_util.safe_set(TestID,(3,2,),945)
# temp/beyler_hood.m:44
TestID = smop_util.safe_set(TestID,(3,3,),947)
# temp/beyler_hood.m:45
TestID = smop_util.safe_set(TestID,(3,4,),951)
# temp/beyler_hood.m:46
Fuel = smop_util.safe_set(Fuel,(4,),'propane')
# temp/beyler_hood.m:48
TestID = smop_util.safe_set(TestID,(4,0,),232)
# temp/beyler_hood.m:50
TestID = smop_util.safe_set(TestID,(4,1,),257)
# temp/beyler_hood.m:51
TestID = smop_util.safe_set(TestID,(4,2,),287)
# temp/beyler_hood.m:52
TestID = smop_util.safe_set(TestID,(4,3,),303)
# temp/beyler_hood.m:53
TestID = smop_util.safe_set(TestID,(4,4,),307)
# temp/beyler_hood.m:54
TestID = smop_util.safe_set(TestID,(4,5,),318)
# temp/beyler_hood.m:55
TestID = smop_util.safe_set(TestID,(4,6,),322)
# temp/beyler_hood.m:56
TestID = smop_util.safe_set(TestID,(4,7,),334)
# temp/beyler_hood.m:57
TestID = smop_util.safe_set(TestID,(4,8,),355)
# temp/beyler_hood.m:58
TestID = smop_util.safe_set(TestID,(4,9,),359)
# temp/beyler_hood.m:59
TestID = smop_util.safe_set(TestID,(4,10,),371)
# temp/beyler_hood.m:60
TestID = smop_util.safe_set(TestID,(4,11,),389)
# temp/beyler_hood.m:61
TestID = smop_util.safe_set(TestID,(4,12,),429)
# temp/beyler_hood.m:62
TestID = smop_util.safe_set(TestID,(4,13,),433)
# temp/beyler_hood.m:63
TestID = smop_util.safe_set(TestID,(4,14,),445)
# temp/beyler_hood.m:64
Fuel = smop_util.safe_set(Fuel,(5,),'propylene')
# temp/beyler_hood.m:66
TestID = smop_util.safe_set(TestID,(5,0,),780)
# temp/beyler_hood.m:68
TestID = smop_util.safe_set(TestID,(5,1,),805)
# temp/beyler_hood.m:69
TestID = smop_util.safe_set(TestID,(5,2,),859)
# temp/beyler_hood.m:70
TestID = smop_util.safe_set(TestID,(5,3,),870)
# temp/beyler_hood.m:71
TestID = smop_util.safe_set(TestID,(5,4,),882)
# temp/beyler_hood.m:72
TestID = smop_util.safe_set(TestID,(5,5,),886)
# temp/beyler_hood.m:73
TestID = smop_util.safe_set(TestID,(5,6,),910)
# temp/beyler_hood.m:74
Fuel = smop_util.safe_set(Fuel,(6,),'toluene')
# temp/beyler_hood.m:76
TestID = smop_util.safe_set(TestID,(6,0,),160)
# temp/beyler_hood.m:78
TestID = smop_util.safe_set(TestID,(6,1,),162)
# temp/beyler_hood.m:79
TestID = smop_util.safe_set(TestID,(6,2,),165)
# temp/beyler_hood.m:80
TestID = smop_util.safe_set(TestID,(6,3,),166)
# temp/beyler_hood.m:81
TestID = smop_util.safe_set(TestID,(6,4,),170)
# temp/beyler_hood.m:82
Species=np.asarray([], dtype='object')
# temp/beyler_hood.m:84
Species = smop_util.safe_set(Species,(0,),'O$_{2}$')
# temp/beyler_hood.m:85
Species = smop_util.safe_set(Species,(1,),'CO$_{2}$')
# temp/beyler_hood.m:86
Species = smop_util.safe_set(Species,(2,),'H$_{2}$O')
# temp/beyler_hood.m:87
Species = smop_util.safe_set(Species,(3,),'CO')
# temp/beyler_hood.m:88
Species = smop_util.safe_set(Species,(4,),'UHC')
# temp/beyler_hood.m:89
Species = smop_util.safe_set(Species,(5,),'Soot')
# temp/beyler_hood.m:90
SaveName=np.asarray([], dtype='object')
# temp/beyler_hood.m:92
SaveName = smop_util.safe_set(SaveName,(0,),'O2')
# temp/beyler_hood.m:93
SaveName = smop_util.safe_set(SaveName,(1,),'CO2')
# temp/beyler_hood.m:94
SaveName = smop_util.safe_set(SaveName,(2,),'H2O')
# temp/beyler_hood.m:95
SaveName = smop_util.safe_set(SaveName,(3,),'CO')
# temp/beyler_hood.m:96
SaveName = smop_util.safe_set(SaveName,(4,),'UHC')
# temp/beyler_hood.m:97
SaveName = smop_util.safe_set(SaveName,(5,),'Soot')
# temp/beyler_hood.m:98
NumPoints=np.asarray([], dtype='object')
# temp/beyler_hood.m:100
NumPoints = smop_util.safe_set(NumPoints,(0,),5)
# temp/beyler_hood.m:101
NumPoints = smop_util.safe_set(NumPoints,(1,),5)
# temp/beyler_hood.m:102
NumPoints = smop_util.safe_set(NumPoints,(2,),5)
# temp/beyler_hood.m:103
NumPoints = smop_util.safe_set(NumPoints,(3,),5)
# temp/beyler_hood.m:104
NumPoints = smop_util.safe_set(NumPoints,(4,),15)
# temp/beyler_hood.m:105
NumPoints = smop_util.safe_set(NumPoints,(5,),7)
# temp/beyler_hood.m:106
NumPoints = smop_util.safe_set(NumPoints,(6,),5)
# temp/beyler_hood.m:107
N_Fuels=7
# temp/beyler_hood.m:109
N_Species=6
# temp/beyler_hood.m:110
X_leg_pos=np.asarray([0.55,0.3,0.2,0.2], dtype='object')
# temp/beyler_hood.m:112
Y_leg_pos=np.asarray([0.55,0.3,0.2,0.2], dtype='object')
# temp/beyler_hood.m:113
# Color per fuel
color=np.asarray([], dtype='object')
# temp/beyler_hood.m:116
color = smop_util.safe_set(color,(0,),'k')
# temp/beyler_hood.m:117
color = smop_util.safe_set(color,(1,),'r')
# temp/beyler_hood.m:118
color = smop_util.safe_set(color,(2,),'b')
# temp/beyler_hood.m:119
color = smop_util.safe_set(color,(3,),'g')
# temp/beyler_hood.m:120
color = smop_util.safe_set(color,(4,),'m')
# temp/beyler_hood.m:121
color = smop_util.safe_set(color,(5,),'c')
# temp/beyler_hood.m:122
color = smop_util.safe_set(color,(6,),'k')
# temp/beyler_hood.m:123
# Marker per fuel
marker=np.asarray([], dtype='object')
# temp/beyler_hood.m:126
marker = smop_util.safe_set(marker,(0,),'o')
# temp/beyler_hood.m:127
marker = smop_util.safe_set(marker,(1,),'s')
# temp/beyler_hood.m:128
marker = smop_util.safe_set(marker,(2,),'d')
# temp/beyler_hood.m:129
marker = smop_util.safe_set(marker,(3,),'^')
# temp/beyler_hood.m:130
marker = smop_util.safe_set(marker,(4,),'v')
# temp/beyler_hood.m:131
marker = smop_util.safe_set(marker,(5,),'>')
# temp/beyler_hood.m:132
marker = smop_util.safe_set(marker,(6,),'<')
# temp/beyler_hood.m:133
from plot_style import *
Marker_Size=7
# temp/beyler_hood.m:136
# Collect data

ExpPlot=np.asarray([], dtype='object')
# temp/beyler_hood.m:140
FDSPlot=np.asarray([], dtype='object')
# temp/beyler_hood.m:141
for f in range(1,N_Fuels+1):
    for s in range(1,NumPoints[f-1]+1):
        FDS_File=outdir + 'Beyler_Hood_' + Fuel[f-1] + '_' + str(TestID[f-1,s-1]) + '_lr_devc.csv'
# temp/beyler_hood.m:144
        fds_data=np.nan_to_num(np.genfromtxt(FDS_File, delimiter = ',')[2:,:], copy = False)
# temp/beyler_hood.m:145
        n_fds=fds_data.shape[1-1]
# temp/beyler_hood.m:146
        for ns in range(1,N_Species+1):
            ExpPlot = smop_util.safe_set(ExpPlot,(f-1,s-1,ns-1,),exp_data[s-1,np.dot((f - 1),N_Species) + ns-1])#NOTE THERE IS NO SECOND -1 FOR THE f
# temp/beyler_hood.m:148
            FDSPlot = smop_util.safe_set(FDSPlot,(f-1,s-1,ns-1,),np.mean(fds_data[n_fds - 60-1:n_fds,1 + ns-1]))
# temp/beyler_hood.m:149

Xmax=np.asarray([0.25,0.25,0.2,0.1,0.1,0.02], dtype='object')
# temp/beyler_hood.m:153
hf=np.asarray([], dtype='object')
# temp/beyler_hood.m:154
XLegendStr=np.asarray([], dtype='object')
# temp/beyler_hood.m:155
hX=np.asarray([], dtype='object')
# temp/beyler_hood.m:156
for ns in range(1,N_Species+1):
    #set(gca,'FontName',Font_Name)
    plt.rcParams["font.family"] = Font_Name
    hf = smop_util.safe_set(hf,(ns-1,),plt.figure(ns))
# temp/beyler_hood.m:158
    #Xmax = max(max(FDSPlot(:,:,ns)));
   #Xmax = max(max(max(ExpPlot(:,:,ns))),Xmax);
   #Xmax = ceil(Xmax*10)/10;
    for f in range(1,N_Fuels+1):
        for s in range(1,NumPoints[f-1]+1):
            XLegendStr = smop_util.safe_set(XLegendStr,(f-1,),Fuel[f-1])
# temp/beyler_hood.m:165
            #hX = smop_util.safe_set(hX,(f-1,),plt.plot(ExpPlot[f-1,s-1,ns-1],FDSPlot[f-1,s-1,ns-1]))
            hX = smop_util.safe_set(hX,(f-1,),plt.plot(ExpPlot[f-1,s-1,ns-1],FDSPlot[f-1,s-1,ns-1],marker=marker[f-1],markersize=Marker_Size,markeredgecolor=color[f-1],markerfacecolor='none',linewidth=Line_Width,linestyle='None')[0])
# temp/beyler_hood.m:167
            #set(hX[f-1],'Marker',marker[f-1],'MarkerSize',Marker_Size,'MarkerEdgeColor',color[f-1],'MarkerFaceColor','none','LineWidth',Line_Width,'LineStyle','none')
            #hold('on')
    xmin=0.0
# temp/beyler_hood.m:178
    ymin=0.0
# temp/beyler_hood.m:179
    xmax=Xmax[ns-1]
# temp/beyler_hood.m:180
    ymax=xmax
# temp/beyler_hood.m:181
    plt.plot(np.asarray([xmin,xmax], dtype='object'),np.asarray([ymin,ymax], dtype='object'),linewidth=.1)
    #axis(np.asarray([xmin,xmax,ymin,ymax], dtype='object'))
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    #set(gca,'PlotBoxAspectRatio',np.asarray([1,1,1], dtype='object'))
    plt.gca().set_aspect('equal')
    xtitle='Measured ' + Species[ns-1] + ' (Mass Fraction)'
# temp/beyler_hood.m:187
    ytitle='Predicted ' + Species[ns-1] + ' (Mass Fraction)'
# temp/beyler_hood.m:188
    #xlabel(xtitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
    plt.xlabel(xtitle,fontsize=Scat_Label_Font_Size)
    #ylabel(ytitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
    plt.ylabel(ytitle,fontsize=Scat_Label_Font_Size)
    #lh=legend(hX,XLegendStr,'Location','NorthWest')
    plt.legend(hX,XLegendStr,loc='upper left',fontsize=Key_Font_Size)
# temp/beyler_hood.m:191
    #set(lh,'FontSize',Key_Font_Size)
    # add VerStr if file is available
    git_file=outdir + 'Beyler_Hood_acetone_117_lr_git.txt'
# temp/beyler_hood.m:196
    from addverstr import *
    addverstr(plt.gca(),git_file,'linear')
    # print to pdf
    #set(gcf,'Visible',Figure_Visibility)
    #set(gcf,'Units',Paper_Units)
    #set(gcf,'PaperUnits',Paper_Units)
    #set(gcf,'PaperSize',np.asarray([Scat_Paper_Width,Scat_Paper_Height], dtype='object'))
    plt.gcf().set_size_inches(Scat_Paper_Width,Scat_Paper_Height)
    #set(gcf,'Position',np.asarray([0,0,Scat_Paper_Width,Scat_Paper_Height], dtype='object'))
    plotname='Beyler_Hood_' + SaveName[ns-1] + '.pdf'
# temp/beyler_hood.m:205
    #print_(gcf,'-dpdf',plotname)
    plt.savefig(plotname,format='pdf')
#plt.show()
