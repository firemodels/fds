import numpy as np
import os.path
import matplotlib.pyplot as plt
# McDermott
# 6-23-2015
# addverstr.m


def addverstr(ha=None,fn=None,pt=None,varargin=None):

    VerStr_Scale_X=0.6
# temp/addverstr.m:7
    VerStr_Scale_Y=1.05
# temp/addverstr.m:8
    Font_Name='times'
# temp/addverstr.m:9
    Font_Interpreter='TeX'
# temp/addverstr.m:10
    try:
        nvararg=varargin.size
# temp/addverstr.m:12
        if nvararg >= 1:
            VerStr_Scale_X=varargin[0]
# temp/addverstr.m:14
        
        if nvararg >= 2:
            VerStr_Scale_Y=varargin[1]
# temp/addverstr.m:15
        
        if nvararg >= 3:
            Font_Name=varargin[2]
# temp/addverstr.m:16
        
        if nvararg >= 4:
            Font_Interpreter=varargin[3]
# temp/addverstr.m:17
    except:
        pass
    
    if os.path.isfile(fn):
        with open(fn,'r') as f:
            VerStr = f.read().strip()
# temp/addverstr.m:20
        x_lim=ha.get_xlim()
# temp/addverstr.m:21
        y_lim=ha.get_ylim()
# temp/addverstr.m:22
        if pt == 'loglog':
            X_VerStr_Position=10 ** (log10(x_lim[0]) + np.dot(VerStr_Scale_X,(log10(x_lim[1]) - log10(x_lim[0]))))
# temp/addverstr.m:24
            Y_VerStr_Position=10 ** (log10(y_lim[0]) + np.dot(VerStr_Scale_Y,(log10(y_lim[1]) - log10(y_lim[0]))))
# temp/addverstr.m:25
        else:
            if pt == 'semilogx':
                X_VerStr_Position=10 ** (log10(x_lim[0]) + np.dot(VerStr_Scale_X,(log10(x_lim[1]) - log10(x_lim[0]))))
# temp/addverstr.m:27
                Y_VerStr_Position=y_lim[0] + np.dot(VerStr_Scale_Y,(y_lim[1] - y_lim[0]))
# temp/addverstr.m:28
            else:
                if pt == 'semilogy':
                    X_VerStr_Position=x_lim[0] + np.dot(VerStr_Scale_X,(x_lim[1] - x_lim[0]))
# temp/addverstr.m:30
                    Y_VerStr_Position=10 ** (log10(y_lim[0]) + np.dot(VerStr_Scale_Y,(log10(y_lim[1]) - log10(y_lim[0]))))
# temp/addverstr.m:31
                else:
                    X_VerStr_Position=x_lim[0] + np.dot(VerStr_Scale_X,(x_lim[1] - x_lim[0]))
# temp/addverstr.m:33
                    Y_VerStr_Position=y_lim[0] + np.dot(VerStr_Scale_Y,(y_lim[1] - y_lim[0]))
# temp/addverstr.m:34
        plt.text(X_VerStr_Position,Y_VerStr_Position,VerStr, fontsize=10, fontname=Font_Name, usetex=True)
        #if isnumeric(VerStr):
            #text(X_VerStr_Position,Y_VerStr_Position,'VerStr ' + str(VerStr),'FontSize',10,'FontName',Font_Name,'Interpreter','none')
        #else:
            #if ischar(VerStr[0]):
                #text(X_VerStr_Position,Y_VerStr_Position,np.asarray([VerStr], dtype='object'),'FontSize',10,'FontName',Font_Name,'Interpreter','none')