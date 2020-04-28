#!/usr/bin/python
#McDermott
#8-30-2016

import subprocess

# read list of directories to move
f = open("dirs.txt")
folders = f.read().splitlines()
#print folders

for i in range(len(folders)):
    #print folders[i]
    #subprocess.call(['ls','-l','Validation'])
    subprocess.call(['git mv '+folders[i]+'/Experimental_Data/* '+folders[i]+'/.'], shell=True)












