#!/usr/bin/python3
#McDermott
#2019-10-08

import struct
import numpy as np

def slread( fname, Tstart, Tend, *args ):

    """
    Reads FDS slice file
    Based on slread.m by Simo Hostikka
    https://github.com/firemodels/fds/blob/master/Utilities/Matlab/scripts/slread.m

    (QQ,Time)=slread(fname,Tstart,Tend [,Tstep]);
      Tstart  is start time
      Tend    is end time
      Nframes is number of slice file frames (FDS default is 1000 between Tstart and Tend)

      QQ      contains the data
      Time    contains the time points
    """

    print(fname)
    print(Tstart)
    print(Tend)

    if len(args)==0:
        Nframes = 1000
    else:
        Nframes = args[0]

    f = open(fname,'rb')

    f.read(4)
    Str1 = f.read(30)
    print(Str1)
    f.read(8)
    Str2 = f.read(30)
    print(Str2)
    f.read(8)
    Str3 = f.read(30)
    print(Str3)
    f.read(8)
    Indx = struct.unpack('6i',f.read(24))
    print(Indx)
    f.read(4)

    # allocate arrays for QQ and Time

    Isize = Indx[1]-Indx[0]+1
    Jsize = Indx[3]-Indx[2]+1
    Ksize = Indx[5]-Indx[4]+1
    if Isize==1:
       M = Jsize
       N = Ksize
    elif Jsize==1:
       M = Isize
       N = Ksize
    else:
       M = Isize
       N = Jsize

    Nframes = max(1,Nframes)
    QQ = np.zeros((N,M,Nframes+1))
    Time = np.zeros(Nframes+1)

    st = 0

    while Time[st] < Tstart:

        f.read(4)
        Time_list = struct.unpack('f',f.read(4))
        Time[st] = Time_list[0]
        f.read(8)
        for m in range(M):
            for n in range(N):
                QQ_list = struct.unpack('f',f.read(4))
                QQ[m,n,st] = QQ_list[0]
        f.read(4)

    while Time[st] < Tend:

        f.read(4)
        Time_list = struct.unpack('f',f.read(4))
        Time[st] = Time_list[0]
        f.read(8)
        for m in range(M):
            for n in range(N):
                QQ_list = struct.unpack('f',f.read(4))
                QQ[m,n,st] = QQ_list[0]
        f.read(4)
        st = st + 1
        if st>Nframes:
            break


    return(QQ,Time)



