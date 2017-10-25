#me -*- coding: utf-8 -*-

#####################################################################################
# Initial import
#####################################################################################
import sys
from pylab import *
from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib.font_manager as fnt


def read_csv(chid, time, quan, nquan, title, isim):

    time_sim  = []
    quan_sim  = []

    for iq in range(nquan):
       quan_sim.append([])

    chid_sim = chid[isim]
    name_chid= "cases/%s/%s_devc.csv" % (chid_sim, chid_sim)
    print 'reading %s'%name_chid

    chid_in  = open(name_chid ,'r')
    chid_input  = chid_in.readlines()
    chid_in.close()

    num=0
    start = 2

    for line in chid_input:

        if start == 1:
           line, null = line.split ("\n")
           titles = line.split (",")
        
        if (start > 0):
            start -= 1
            continue
        num=num+1

        line, null = line.split ("\n")
        values = line.split (",")

        t = float(values[0])
        time_sim.append(t)

        for iq in range(nquan):
           q = float(values[iq+1])
           quan_sim[iq].append(q)
   
    time.append(time_sim)
    quan.append(quan_sim)
    title.append(titles)


def plot_csv(chid, time, quan, iquan, title, nsim, obsttype):

    name = '%s_%s'%(obsttype, title[0][iquan+1].strip('\"'))

    fig = plt.figure (facecolor='w')
    ax = fig.add_subplot(111)

    markers = ["s","D","o","v","h","^","<","s","D","o","v","h","^","<"]
    colors  = ["c","m","y","k","r","b","g","r","m","y","k","c","b","g"]

    legsize = fnt.FontProperties(size=8)

    index  = []
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

    pstart0 = max(pstart[:])
    pend0   = min(pend[:])

    #print pstart
    #print pend
    #print pstart0
    #print pend0

    step = (pend0 - pstart0)/(nsim+1)
 
    min_val=[]
    max_val=[]

    for isim in range(nsim):
       #print i, min(quan[i]), max(quan[i])
       #print 'isim=',isim, pstart0, pend0, len(quan[isim][iquan])
       min_val.append(min(quan[isim][iquan][pstart0:pend0]))
       max_val.append(max(quan[isim][iquan][pstart0:pend0]))
       #min_val.append(min(quan[isim][iquan]))
       #max_val.append(max(quan[isim][iquan]))

    min_total = min(min_val)
    max_total = max(max_val)
    #print 'min_total:',min_total
    #print 'max_total:',max_total

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85])

    for isim in range(nsim):
      pstart0 += step
      mark =  (pstart0, pend0)
      ax.plot (time[isim], quan[isim][iquan],'-r', linewidth=1.0, markevery=mark, marker=markers[isim], color = colors[isim])
      legend.append(chid[isim])

    ax.legend(legend, prop=legsize, loc="lower center", 
              bbox_to_anchor=(0.50, -0.25), fancybox=True, shadow=True, ncol=3)

    ax.grid(True)
    #ax.set_title('%s meshes'%(top),fontsize=30)

    #ax.set_xlabel('Time [s]',fontsize=16)
    #ax.set_ylabel(r'%s' %title[0][iquan+1].strip('\"'),fontsize=16)


    ax.set_xlabel('Time [s]',fontsize=16)
    yname='none'
    if  'PRES_ITE' in name:
       yname = 'Pressure iterations'
    elif 'PRES' in name:
       yname = 'Pressure'
    elif  'UVEL' in name:
       yname = 'U-velocity'
    elif  'VVEL' in name:
       yname = 'V-velocity'
    elif  'WVEL' in name:
       yname = 'W-velocity'
    elif  'SCARC_ITE' in name:
       yname = 'ScaRC iterations'
    elif  'SCARC_CONV' in name:
       yname = 'ScaRC convergence rate'
    elif  'SCARC_RES' in name:
       yname = 'ScaRC residual'
    elif  'VEL_ERR' in name:
       yname = 'Velocity error'

    ax.set_ylabel(r'%s' %yname,fontsize=16)



    #ax.set_xlim(tstart,1.1*tend)
    #ax.set_ylim(min_total - abs(min_total*0.1), max_total + abs(max_total*0.06))
    #ax.set_ylim(-0.05, 0.2)

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

    output = "pictures/%s.png" % (name)
    savefig(output)
    #show()

tstart= 0.05
tend  = 1.0


obsttype = sys.argv[1]

nquan = 8
chid=[]


discretizations = ['structured']
#topologies      = [[1,1,1],[2,1,2]]
topologies      = [[1,1,1]]
cells           = [[16,1,16]]
solvers         = ['fft','cg_ssor','cg_fft','gmg_ssor','fft_cg_ssor']
tolerances      = [6]


for discret in discretizations:
   for top in topologies[:]:
      nmeshes = int(top[0])*int(top[1])*int(top[2])
      for solver in solvers[:]:
         if discret == 'structured':
            for tol in tolerances:
               chid.append('s%d_%s_tol%d_%s' %(nmeshes, obsttype, tol, solver))
         elif discret == 'unstructured':
            if 'fft' in solver: continue
            chid.append('u%d_%s_tol%d_%s' %(nmeshes, obsttype, tol, solver))

nsim = len(chid)

time = []
quan = []
title = []

for isim in range(nsim):
   read_csv(chid, time, quan, nquan, title, isim)

for iquan in range(nquan):
   plot_csv(chid, time, quan, iquan, title, nsim, obsttype)





