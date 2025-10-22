#McDermott
#2019-10-08

import struct
import numpy as np

def slread( fname, Tstart, Tend, *args, **kwargs ):

    """
    Reads FDS slice file
    Based on slread.m by Simo Hostikka
    https://github.com/firemodels/fds/blob/master/Utilities/Matlab/scripts/slread.m
    (QQ,Time)=slread(fname,Tstart,Tend [,Nframes, gridskip, timeskip]);
      Tstart  is start time
      Tend    is end time
      Nframes is number of slice file frames (FDS default is 1000 between Tstart and Tend)
      gridskip is a skip rate for reading cells: =2 --> read every other cell, etc.
      timeskip is a skip rate for frames: =2 --> read every other frame, etc.
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
        
    gridskip = kwargs['gridskip'] if 'gridskip' in kwargs else 1
    timeskip = kwargs['timeskip'] if 'timeskip' in kwargs else 1
        
    f = open(fname,'rb')

    f.read(4)
    Str1 = f.read(30)       # quantity
    print(Str1)
    f.read(8)
    Str2 = f.read(30)       # short name
    print(Str2)
    f.read(8)
    Str3 = f.read(30)       # units
    print(Str3)
    f.read(8)
    Indx = struct.unpack('6i',f.read(24))  # index bounds i,j,k
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
    QQ = np.zeros((M,N))
    Time = np.zeros(Nframes+1)
    
    ii = np.arange(0,M,gridskip)
    jj = np.arange(0,N,gridskip)
    tt = np.arange(0,Nframes+1,timeskip)
    Qskip  = np.zeros((len(ii),len(jj),len(tt)))

    st = 0

    while Time[st] < Tstart:

        f.read(4)
        Time_list = struct.unpack('f',f.read(4))
        Time[st] = Time_list[0]
        f.read(8)
        for n in range(N):
            for m in range(M):
                QQ_list = struct.unpack('f',f.read(4))
                QQ[m,n] = QQ_list[0]
        f.read(4)

    while Time[st] < Tend:

        f.read(4)
        Time_list = struct.unpack('f',f.read(4))
        Time[st] = Time_list[0]
        f.read(8)
        for n in range(N):
            for m in range(M):
                QQ_list = struct.unpack('f',f.read(4))
                QQ[m,n] = QQ_list[0]
        if st%timeskip==0:
            print(st,int(st/timeskip))
            Qskip[:,:,int(st/timeskip)] = QQ[np.ix_(ii,jj)]
        f.read(4)
        st = st + 1
        if st>Nframes:
            break


    return(Qskip,Time[tt])

