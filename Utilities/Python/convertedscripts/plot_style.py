import numpy as np
import matplotlib.pyplot as plt
# K. McGrattan
# 7-16-2009
# plot_style.m

# Preferred style for FDS and Smokeview plots

# Usage: When creating a new figure use the following:

# plot_style
# figure
# set(gca,'Units',Plot_Units)
# set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

# ==>> plot(X,Y,etc)
# ==>> xlabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
# ==>> ylabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

# set(gca,'FontName',Font_Name)
# set(gca,'FontSize',Label_Font_Size)

# ==>> lh=legend(...);
# ==>> set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

# set(gcf,'Visible',Figure_Visibility);
# set(gcf,'Units',Paper_Units);
# set(gcf,'PaperUnits',Paper_Units);
# set(gcf,'PaperSize',[Paper_Width Paper_Height]);
# set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
# print(gcf,'-dpdf',plotname);

style='fds'
# temp/plot_style.m:31

if 'fds' == style:
    # Font properties
    Font_Name='times'
# temp/plot_style.m:38
    Font_Interpreter='TeX'
    plt.rcParams["text.usetex"] = True
# temp/plot_style.m:39
    Key_Font_Size=11
# temp/plot_style.m:40
    Title_Font_Size=16
# temp/plot_style.m:41
    Label_Font_Size=16
# temp/plot_style.m:42
    Scat_Title_Font_Size=14
# temp/plot_style.m:43
    Scat_Label_Font_Size=14
# temp/plot_style.m:44
    Marker_Size=4
# temp/plot_style.m:45
    Line_Width=1.0
# temp/plot_style.m:48
    Plot_Units='inches'
# temp/plot_style.m:51
    Plot_Width=5.0
# temp/plot_style.m:52
    Plot_Height=3.4
# temp/plot_style.m:53
    Plot_X=1.2
# temp/plot_style.m:54
    Plot_Y=0.8
# temp/plot_style.m:55
    Scat_Plot_Width=4.75
# temp/plot_style.m:57
    Scat_Plot_Height=4.75
# temp/plot_style.m:58
    Scat_Plot_X=1.0
# temp/plot_style.m:59
    Scat_Plot_Y=0.75
# temp/plot_style.m:60
    Subtitle_Text_Offset=0.05
# temp/plot_style.m:61
    VerStr_Scale_X=0.6
# temp/plot_style.m:63
    VerStr_Scale_Y=1.05
# temp/plot_style.m:64
    Paper_Units='inches'
# temp/plot_style.m:67
    Paper_Width=6.5
# temp/plot_style.m:68
    Paper_Height=4.5
# temp/plot_style.m:69
    Scat_Paper_Height=6.0
# temp/plot_style.m:70
    Scat_Paper_Width=6.0
# temp/plot_style.m:71
    Figure_Visibility='off'
# temp/plot_style.m:74
    Image_File_Type='-dpdf'
# temp/plot_style.m:75
else:
    if 'paper' == style:
        # Font properties
        Font_Name='helvetica'
# temp/plot_style.m:80
        Font_Interpreter='TeX'
        plt.rcParams["text.usetex"] = True
# temp/plot_style.m:81
        Key_Font_Size=12
# temp/plot_style.m:82
        Title_Font_Size=20
# temp/plot_style.m:83
        Label_Font_Size=20
# temp/plot_style.m:84
        Marker_Size=10
# temp/plot_style.m:85
        Scat_Title_Font_Size=14
# temp/plot_style.m:86
        Scat_Label_Font_Size=14
# temp/plot_style.m:87
        Line_Width=1
# temp/plot_style.m:90
        # Plot properties
        Plot_Units='normalized'
# temp/plot_style.m:93
        Pos=np.asarray([0.15,0.15,0.775,0.815], dtype='object')
# temp/plot_style.m:94
        YLabel_Offset=- 0.05
# temp/plot_style.m:95
        XLabel_Offset=- 0.05
# temp/plot_style.m:96
        Plot_X=Pos[0]
# temp/plot_style.m:97
        Plot_Y=Pos[1]
# temp/plot_style.m:98
        Plot_Width=Pos[2]
# temp/plot_style.m:99
        Plot_Height=Pos[3]
# temp/plot_style.m:100
        Scat_Plot_Width=4.75
# temp/plot_style.m:102
        Scat_Plot_Height=4.75
# temp/plot_style.m:103
        Scat_Plot_X=0.75
# temp/plot_style.m:104
        Scat_Plot_Y=0.75
# temp/plot_style.m:105
        Subtitle_Text_Offset=0.05
# temp/plot_style.m:106
        VerStr_Scale_X=0.6
# temp/plot_style.m:108
        VerStr_Scale_Y=1.05
# temp/plot_style.m:109
        Paper_Units='inches'
# temp/plot_style.m:112
        Paper_Pos=np.asarray([0.25,2.5,8.0,6.0], dtype='object')
# temp/plot_style.m:113
        Paper_Width=Paper_Pos[2]
# temp/plot_style.m:114
        Paper_Height=Paper_Pos[3]
# temp/plot_style.m:115
        Scat_Paper_Height=6.0
# temp/plot_style.m:116
        Scat_Paper_Width=6.0
# temp/plot_style.m:117
        Figure_Visibility='on'
# temp/plot_style.m:120
        Image_File_Type='-dpdf'
# temp/plot_style.m:121