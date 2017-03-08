import os, sys

path =  os.getcwd()
filenames = os.listdir(path)

old = sys.argv[1]
new = sys.argv[2]

for filename in filenames:
    if ".fds" in filename:
        os.rename(filename, filename.replace(old, new))
