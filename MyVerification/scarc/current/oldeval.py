#me -*- coding: utf-8 -*-

#####################################################################################
# Initial import
#####################################################################################
import sys
from pylab import *
from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib.font_manager as fnt
import matplotlib


site = []
sres = []
scon = []

def read_csv(chid, time, pres, uvel, wvel, verr, pite, site, sres, scon, i):

    time0  =[]
    pres0  =[]
    uvel0  =[]
    wvel0  =[]
    verr0  =[]
    pite0  =[]
    site0  =[]
    sres0  =[]
    scon0  =[]

    chid0 = chid[i]
    name_chid= "%s_devc.csv" % (chid0)
    print 'reading %s'%name_chid

    chid_in  = open(name_chid ,'r')
    chid_input  = chid_in.readlines()
    chid_in.close()

    num=0
    start = 3

    for line in chid_input:

        if (start > 0):
            start -= 1
            continue
        num=num+1

        line, null = line.split ("\n")
        quan = line.split (",")

        t = float(quan[0])
        time0.append(t)

        p = float(quan[1])
        pres0.append(p)

        u = float(quan[2])
        uvel0.append(u)

        w = float(quan[3])
        wvel0.append(w)

        ve = float(quan[4])
        verr0.append(ve)

        pi = float(quan[5])
        pite0.append(pi)

        si = float(quan[6])
        site0.append(si)

        sr = float(quan[7])
        sres0.append(sr)

        sc = float(quan[8])
        scon0.append(sc)

   
    time.append(time0)
    pres.append(pres0)
    uvel.append(uvel0)
    wvel.append(wvel0)
    verr.append(verr0)
    pite.append(pite0)
    site.append(site0)
    sres.append(sres0)
    scon.append(scon0)


def plot_csv(chid, time, quan, name, tstart, tend, nmeshes, obst):


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
          if time[i][j]>=tend-0.02:
             pend.append(j)
             break


    nplots = nsim
        
    pstart0   = max(pstart[:])
    pend0     = min(pend[:])

    #print pstart
    #print pend
    #print pstart0
    #print pend0

    step = (pend0 - pstart0)/(nplots+1)
 
    min_val=[]
    max_val=[]

    for i in range(nsim):
       #print 'name[0]=',name[0], chid[i]
       if name[0] == "s" and "cg" not in chid[i]: continue
       #print i, min(quan[i]), max(quan[i])
       #print 'i=',i, pstart0, pend0, len(quan[i])
       min_val.append(min(quan[i][pstart0:pend0]))
       max_val.append(max(quan[i][pstart0:pend0]))
       #min_val.append(min(quan[i]))
       #max_val.append(max(quan[i]))

    min_total = min(min_val)
    max_total = max(max_val)
    #print 'min:',min_val
    #print 'max:',max_val

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85])

    for i in range(nsim):
       if name[0] == "s" and "cg" not in chid[i]: continue
       pstart0 += step
       mark =  (pstart0, pend0)
       ax.plot (time[i], quan[i],'-r', linewidth=1.0, markevery=mark, marker=markers[i], color = colors[i])
       legend.append(chid[i])

    ax.legend(legend, prop=legsize, loc="lower center", 
              bbox_to_anchor=(0.50, -0.25), fancybox=True, shadow=True, ncol=4)

    ax.grid(True)
    #ax.set_title('%s meshes'%(nmeshes),fontsize=30)

    ax.set_xlabel('Time [s]',fontsize=16)
    if name == 'pres':
       yname = 'pressure'
    elif name == 'uvel':
       yname = 'u-velocity'
    elif name == 'vvel':
       yname = 'v-velocity'
    elif name == 'wvel':
       yname = 'w-velocity'
    elif name == 'site':
       yname = 'ScaRC iterations'
    elif name == 'scon':
       yname = 'ScaRC convergence rate'
    elif name == 'sres':
       yname = 'ScaRC residual'
    elif name == 'pite':
       yname = 'pressure iterations'
    elif name == 'verr':
       yname = 'velocity error'
   
    ax.set_ylabel(r'%s' %yname,fontsize=16)


    total = max_total - min_total
    if total == 0.0:
       total_eps = 1E-1
    else:
       total_eps = total*0.1

    ax.set_xlim(tstart,tend)
    ax.set_ylim(min_total - total_eps, max_total + total_eps)
    print 'name=',name,': limits=',min_total, min_total-total_eps, max_total, max_total+total_eps

    labels = ax.get_xticklabels() + ax.get_yticklabels()
    for label in labels:
        label.set_size(13)

    output = "pictures/s%s_%s_%s.png" % (nmeshes, obst, name)
    savefig(output)
    #show()


nmeshes = int(sys.argv[1])
obst    = sys.argv[2]
tstart = 0.1
tend   = 1.00

chid=[]
chid.append('s%d_%s_fft' %(nmeshes, obst))
chid.append('s%d_%s_cg_ssor' %(nmeshes, obst))
chid.append('s%d_%s_cg_fft' %(nmeshes, obst))
chid.append('s%d_%s_gmg_ssor' %(nmeshes, obst))
#chid.append('s%d_%s_fft_cg_ssor' %(nmeshes, obst))


nsim = len(chid)

time = []
pres = []
uvel = []
wvel = []
verr = []
pite = []
site = []
sres = []
scon = []

for i in range(nsim):
   read_csv(chid, time, pres, uvel, wvel, verr, pite, site, sres, scon, i)

plot_csv(chid, time, pres, 'pres', tstart, tend, nmeshes, obst)
plot_csv(chid, time, uvel, 'uvel', tstart, tend, nmeshes, obst)
plot_csv(chid, time, wvel, 'wvel', tstart, tend, nmeshes, obst)
plot_csv(chid, time, verr, 'verr', tstart, tend, nmeshes, obst)
plot_csv(chid, time, pite, 'pite', tstart, tend, nmeshes, obst)
plot_csv(chid, time, site, 'site', tstart, tend, nmeshes, obst)
plot_csv(chid, time, sres, 'sres', tstart, tend, nmeshes, obst)
plot_csv(chid, time, scon, 'scon', tstart, tend, nmeshes, obst)





