#me -*- coding: utf-8 -*-

#####################################################################################
# Initial import
#####################################################################################
import sys
from pylab import *
from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib.font_manager as fnt



def read_csv(chid, time, uvel, vvel, wvel, pres, site, sres, srate, i):

    time0  =[]
    uvel0  =[]
    vvel0  =[]
    wvel0  =[]
    pres0  =[]
    site0  =[]
    sres0  =[]
    srate0 = []

    chid0 = chid[i]
    name_chid= "%s_devc.csv" % (chid0)
    print 'reading %s'%name_chid

    chid_in  = open(name_chid ,'r')
    chid_input  = chid_in.readlines()
    chid_in.close()

    num=0
    start = 2

    for line in chid_input:

        if (start > 0):
            start -= 1
            continue
        num=num+1

        line, null = line.split ("\n")
        quan = line.split (",")

        t = float(quan[0])
        time0.append(t)

        u = float(quan[1])
        uvel0.append(u)

        v = float(quan[2])
        vvel0.append(v)

        w = float(quan[3])
        wvel0.append(w)

        p = float(quan[4])
        pres0.append(p)

        if 'susi' in chid0:

           ite = float(quan[5])
           site0.append(ite)

           res = float(quan[6])
           sres0.append(res)

           rate = float(quan[7])
           srate0.append(rate)

   
    time.append(time0)
    uvel.append(uvel0)
    vvel.append(vvel0)
    wvel.append(wvel0)
    pres.append(pres0)
    site.append(site0)
    sres.append(sres0)
    srate.append(srate0)


def plot_csv(chid, time, quan, name, tstart, tend, nmeshes):


    fig = plt.figure (facecolor='w')
    ax = fig.add_subplot(111)

    markers = ["s","D","o","v","h","^","<"]
    colors  = ["c","m","y","k","r","b","g"]
    colors  = ["r","m","y","k","c","b","g"]

    legsize = fnt.FontProperties(size=8)

    clen = len(chid)
    
    nsim = len(time)

    index = []
    pstart = []
    pend   = []
    legend = []

    for i in range(nsim):
       for j in range(len(time[i])):
          if time[i][j]>=tstart:
             pstart.append(j)
             break

    for i in range(nsim):
       for j in range(len(time[i])):
          if time[i][j]>=tend:
             pend.append(j)
             break


    nplots = nsim
        
    pstart0   = max(pstart[1:])
    pend0     = min(pend[1:])

    print pstart
    print pend
    print pstart0
    print pend0

    step = (pend0 - pstart0)/(nplots+1)
 
    min_val=[]
    max_val=[]

    for i in range(nsim):
       #print i, min(quan[i]), max(quan[i])
       print 'i=',i, pstart0, pend0, len(quan[i])
       #min_val.append(min(quan[i][pstart0:pend0]))
       #max_val.append(max(quan[i][pstart0:pend0]))
       min_val.append(min(quan[i]))
       max_val.append(max(quan[i]))

    min_total = min(min_val)
    max_total = max(max_val)
    min_total = -0.1
    max_total = 1.0
    print 'min:',min_val
    print 'max:',max_val

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85])

    for i in range(nsim):
      pstart0 += step
      mark =  (pstart0, pend0)
      ax.plot (time[i], quan[i],'-r', linewidth=1.0, markevery=mark, marker=markers[i], color = colors[i])
      legend.append(chid[i])

    ax.legend(legend, prop=legsize, loc="lower center", 
              bbox_to_anchor=(0.50, -0.25), fancybox=True, shadow=True, ncol=5)

    ax.grid(True)
    #ax.set_title('%s meshes'%(top),fontsize=30)

    ax.set_xlabel('Time [s]',fontsize=16)
    ax.set_ylabel(r'%s' %name,fontsize=16)

    print min_total, max_total

    ax.set_xlim(tstart,tend)
    ax.set_ylim(min_total - abs(min_total*0.06), max_total + abs(max_total*0.06))
    ax.set_ylim(-0.3, 1.5)

    labels = ax.get_xticklabels() + ax.get_yticklabels()
    for label in labels:
        label.set_size(13)

    output = "pictures/%s_%dmesh_ts%.2f_te%.2f.png" % (name, nmeshes, tstart, tend)
    savefig(output)
    #show()


nmeshes = int(sys.argv[1])
tstart  = float(sys.argv[2])
tend    = float(sys.argv[3])

chid=[]
chid.append('fft_1x1')
chid.append('cg_1x1')
chid.append('mkl_1x1')
chid.append('gmg_1x1')
chid.append('gmg_mkl_1x1')

nsim = len(chid)

time = []
pres = []
uvel = []
vvel = []
wvel = []
site = []
sres = []
srate= []

for i in range(nsim):
   print 'reading ',chid[i]
   read_csv(chid, time, uvel, vvel, wvel, pres, site, sres, srate, i)

plot_csv(chid, time, pres, 'Pressure', tstart, tend, nmeshes)





