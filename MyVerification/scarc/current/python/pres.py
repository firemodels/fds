#me -*- coding: utf-8 -*-

#####################################################################################
# Initial import
#####################################################################################
import sys
from pylab import *
from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib.font_manager as fnt



def read_csv(top, time, pres1, pres2, pres3, site, sres):

    name_top= "%s_devc.csv" % (top)
    print 'reading %s'%name_top

    top_in  = open(name_top ,'r')
    top_input  = top_in.readlines()
    top_in.close()

    num=0
    start = 2

    for line in top_input:

        if (start > 0):
            start -= 1
            continue
        num=num+1

        line, null = line.split ("\n")
        quan = line.split (",")

        t = float(quan[0])
        time.append(t)

        u = float(quan[1])
        pres1.append(u)

        v = float(quan[2])
        pres2.append(v)

        w = float(quan[3])
        pres3.append(w)

        ve = float(quan[4])
        site.append(ve)

        vi = float(quan[5])
        sres.append(vi)

   

def plot_csv3(top, time, pres1, pres2, pres3, name, tstart, tend):


    fig = plt.figure (facecolor='w')
    ax = fig.add_subplot(111)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85])

    ax.plot (time, pres1,'-r', linewidth=1.0)
    ax.plot (time, pres2,'-g', linewidth=1.0)
    ax.plot (time, pres3,'-b', linewidth=1.0)

    legend = ['pres1','pres2','pres3']

    legsize = fnt.FontProperties(size=8)
    ax.legend(legend, prop=legsize, loc="lower center", 
              bbox_to_anchor=(0.50, -0.25), fancybox=True, shadow=True, ncol=3)

    ax.grid(True)
    #ax.set_title('%s meshes'%(top),fontsize=30)

    ax.set_xlabel('Time [s]',fontsize=16)
    ax.set_ylabel(r'%s' %name,fontsize=16)

    ax.set_xlim(tstart,tend)
    #ax.set_ylim(min_total - abs(min_total*0.06), max_total + abs(max_total*0.06))
    #ax.set_ylim(-0.05, 0.2)

    labels = ax.get_xticklabels() + ax.get_yticklabels()
    for label in labels:
        label.set_size(13)

    output = "pictures/%s_%s.png" % (name, top)

    savefig(output)
    #show()


top = sys.argv[1]
tstart=0.0
tend=25.0

print 'reading ',top
time  = []
pres1 = []
pres2 = []
pres3 = []
site  = []
sres  = []
read_csv(top, time, pres1, pres2, pres3, site, sres)

plot_csv3(top, time, pres1, pres2, pres3, 'pres', tstart, tend)





