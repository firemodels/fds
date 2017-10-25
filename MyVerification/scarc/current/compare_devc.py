#me -*- coding: utf-8 -*-

#####################################################################################
# Initial import
#####################################################################################
import sys
from pylab import *
from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib.font_manager as fnt


def read_csv(chid, time, quan, nquan, title):

    time_sim  = []
    quan_sim  = []

    for iq in range(nquan):
       quan_sim.append([])

    chid_sim = chid[0]
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


def plot_csv(chid, time, quan, title):

    name = '%s_all'%(chid[0])

    fig = plt.figure (facecolor='w')
    ax = fig.add_subplot(111)

    markers = ["s","D","o","v","h","^","<","s","D","o","v","h","^","<"]
    colors  = ["c","m","y","k","r","b","g","r","m","y","k","c","b","g"]

    legsize = fnt.FontProperties(size=8)

    index  = []
    pstart = []
    pend   = []
    legend = []

    print "len(time)", len(time[0])
    for j in range(len(time[0])):
       if time[0][j]>=tstart:
          pstart.append(j)
          break

    for j in range(len(time[0])):
       if time[0][j]>=tend:
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

    for iquan in range(nquan):
       min_val.append(min(quan[0][iquan][pstart0:pend0]))
       max_val.append(max(quan[0][iquan][pstart0:pend0]))
       #min_val.append(min(quan[0][iquan]))
       #max_val.append(max(quan[0][iquan]))

    min_total = min(min_val)
    max_total = max(max_val)
    #print 'min:',min_val
    #print 'max:',max_val

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85])

    for iquan in range(nquan):
      pstart0 += step
      mark =  (pstart0, pend0)
      ax.plot (time[0], quan[0][iquan],'-r', linewidth=1.0, markevery=mark, marker=markers[iquan], color = colors[iquan])
      legend.append(title[0][iquan+1].strip('\"'))

    ax.legend(legend, prop=legsize, loc="lower center", 
              bbox_to_anchor=(0.50, -0.25), fancybox=True, shadow=True, ncol=4)

    ax.grid(True)
    #ax.set_title('%s meshes'%(top),fontsize=30)

    ax.set_xlabel('Time [s]',fontsize=16)
    ax.set_ylabel(r'%s' %chid[0],fontsize=16)

    #print min_total, max_total

    ax.set_xlim(tstart,1.1*tend)
    ax.set_ylim(min_total - abs(min_total*0.1), max_total + abs(max_total*0.06))
    #ax.set_ylim(-0.05, 0.2)

    labels = ax.get_xticklabels() + ax.get_yticklabels()
    for label in labels:
        label.set_size(13)

    output = "pictures/%s.png" % (name)
    savefig(output)
    #show()

tstart=0.0
tend  =15.0

chid=[]

zone_type = sys.argv[1]
solver = sys.argv[2]

nquan = 3


discretizations = ['structured']
topologies      = [[1,1,1]]
#cells           = [[11,11,5]]
solvers         = [solver]
tolerances      = [sys.argv[3]]



for discret in discretizations:
   for top in topologies[:]:
      nmeshes = int(top[0])*int(top[1])*int(top[2])
      for solver in solvers[:]:
         if discret == 'structured':
            for tol in tolerances:
               #chid.append('s%d_n%s_%s%s' %(nmeshes, cell, solver, tol))
               chid.append('s%d_%s%s_%s' %(nmeshes, solver, tol, zone_type))
         elif discret == 'unstructured':
            if 'fft' in solver: continue
            chid.append('u%d_%s%s_%s' %(nmeshes, solver, 0, zone_type))

nsim = len(chid)
print 'nsim=',nsim

time = []
quan = []
title = []

read_csv(chid, time, quan, nquan, title)
plot_csv(chid, time, quan, title)





