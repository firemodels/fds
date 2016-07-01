# built by koverholt with an assist from cweinschenk
#!/usr/bin/env python

"""
Display the per-commit size of a file in a git repository.
$ python counter.py setup.py
82e2452, 2016-02-23T23:57:25-05:00, 14 lines
94e26ae, 2016-02-24T12:01:39-06:00, 15 lines
26e51c5, 2016-03-07T13:35:19-06:00, 19 lines
fed0eb8, 2016-03-15T15:03:36-05:00, 18 lines
8b996ee, 2016-03-22T13:45:43-04:00, 18 lines
6c1c996, 2016-03-22T20:04:41-04:00, 53 lines
336f19e, 2016-04-07T15:01:09-05:00, 43 lines
ee9bd08, 2016-04-25T13:30:58-04:00, 41 lines
"""
import subprocess
import re
import sys
import numpy as np
import pandas as pd
from pylab import *
import matplotlib.pyplot as plt

def git_line_count(fn):
    data = []
    git = subprocess.Popen(["git", "log", "--shortstat", "--reverse",
                            "--pretty=format:%h %at", fn],
                           stdout=subprocess.PIPE)
    out, err = git.communicate()
    out = out.decode('ascii')

    line_count = 0
    for line in out.split('\n'):
        if not line.strip():
            continue
        if line[0] != ' ':
            hash, time = line.split(' ')
        elif 'file' in line:
            try:
                insertions = int(re.findall(r'(\d+) insertion[\w]?\(\+\)', line)[0])
                line_count += insertions
            except IndexError:
                pass

            try:
                deletions = int(re.findall(r'(\d+) deletion[\w]?\(\-\)', line)[0])
                line_count -= deletions
            except IndexError:
                pass

            # print("%s, %s, %d lines" % (hash, time, line_count))
            data.append([hash,time,line_count])
    return(data)

file = '../../Verification/FDS_Cases.sh'
gitdata = git_line_count(file)
data = np.asarray(gitdata)
df = pd.DataFrame({'date':data[:,1],'num_lines':data[:,2]},dtype=int)
dates = pd.to_datetime(df.date.values,unit='s')

#convert datetime
conv = np.vectorize(strpdate2num('%Y-%m-%d'))

fig = plt.figure()
plt.plot(dates,.92*df['num_lines'],lw=2)
plt.axvline(conv('2012-06-26'),linestyle='-',color = '#000000')
plt.text(conv('2012-06-30'), 210, 'Firebot Starts', 
     horizontalalignment='left',
     verticalalignment='center')
plt.axvline(conv('2013-11-04'),linestyle='-',color = '#000000')
plt.text(conv('2013-11-08'), 220, 'FDS 6.0.0 Released', 
     horizontalalignment='left',
     verticalalignment='center')
plt.axvline(conv('2016-04-06'),linestyle='-',color = '#000000')
plt.text(conv('2016-04-10'), 470, 'Latest Release: FDS 6.4.0', 
     horizontalalignment='right',
     verticalalignment='center')
xlabel('Time')
ylabel('FDS Verification Cases')
savefig('Verification_Cases.pdf')
