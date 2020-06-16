###############################################################################
# David O. Lignell
# May 21, 2020

# get_var_slice.py function calls and extends the slread.py script. The
# slread.py script reads a given slice file. This will be for a particular
# mesh, for a particular variable. The get_var_slice.py script here provides a
# simple interface that will read all mesh files for a given variable for a
# given slice and combine them into a single variable for easy plotting or
# analysis. This is convenient for computing average or fluctuating fields, or
# otherwise processing the data contained in a slice. The interface requires
# minimal information from the user and reads the FDS input file to gather
# required information.
#
# A complete test case setup is available on Github at
# https://github.com/BYUignite/fds_read_slice.
#
# Note, if this file is called as main, e.g., as: python3 get_var_slice.py,
# then it will call the get_var_slice function and plot the results. See the
# bottom. You'll need to edit key vars below.
# 
# A Jupyter notebook plot_slice.ipynb driver that calls the main function get_var_slice and
# plots the results is provided.
#
# Use Python 3

from slread import *
import numpy as np
import matplotlib.pyplot as plt
import re
from typing import List, Tuple
import sys, os

###############################################################################

def get_var_slice(inputFileName:str, varName:str = 'TEMPERATURE',
                  planeNumber:int=1,
                  gridskip:int=1,
                  timeskip:int=1) -> Tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    '''
    Input: inputFileName (*.fds file)
    Input: varName (should match variable names in slice files)
    Input: planeNumber (which "plane" is the variable in; which occurence is the var in input file)
    Input: gridskip (e.g., =2 to skip every other grid point)
    Input: timeskip (e.g., =2 to skip every other time point)
    Return: variable, with X grid, Y grid for contour plots, and plane direction for reference.
    '''
    
    #---------- open input file, get list of lines
    
    with open(inputFileName) as ifile:    
        iflines = ifile.readlines()
        iflines = _join_continued_lines(iflines)
        
        dir = re.search('/',inputFileName)
        if dir:
            dir = re.search('(\S*)/', inputFileName).group(1)
        else:
            dir = "./"
            
    #---------- get chid, tstart, tend, nframes from input file
    
    tstart = 0
    for line in iflines:
        if line.find('&TIME') == 0:
            tend = float(re.search('T_END\s*=\s*(\d*\.?\d*)', line).group(1))
        if line.find('&DUMP') == 0:
            nframes = int(re.search('NFRAMES\s*=\s*(\d*)', line).group(1))
        if line.find('&HEAD') == 0:
            chid = re.search("CHID\s*=\s*'(.+?)'", line).group(1)
            
    #nframes = 126 # doldb
            
        
    
    #---------- get varprops (list of tuples with: varname, plane dir, plane loc),
    #---------- also get index of variable of interest in this list: varind
    
    varprops = _get_varprops(iflines)
    
    icount = 0
    varind = -1
    for i,var in enumerate(varprops):
        if var[0] == varName:
            icount += 1
            if icount == planeNumber:
                varind = i
                break
                
    #---------- get mesh topology and list of meshes for given var in plane
    
    mx, my, mz, nmx, nmy, nmz = _get_mesh_topology(iflines)
    
    meshes, ivarinmesh = _get_meshes_ivarinmesh(varind, varprops, mx,my,mz)
    
    #---------- read var from various mesh slice files, concatenate into single plane
    
    if varprops[varind][1]=='X':
        nm1 = nmy
        nm2 = nmz
    elif varprops[varind][1]=='Y':
        nm1 = nmx
        nm2 = nmz
    else:       # =='Z'
        nm1 = nmx
        nm2 = nmy
        
    for j in range(nm2):
        for i in range(nm1):
            imesh = j*nm1+i     # index in list of meshes
            fname = f"{dir}/{chid}_{meshes[imesh]+1:04d}_{ivarinmesh[imesh]+1:02d}.sf"
            with _hide_print():
                var, times = slread(fname, tstart, tend, nframes, gridskip=gridskip, timeskip=timeskip)
            if i==0:
                v1 = var.copy()
            else:
                v1 = np.concatenate((v1,var), axis=0)
        if j==0:
            V = v1.copy()
        else:
            V = np.concatenate((V,v1), axis=1)
            
    #-------------
    
    planedir = varprops[varind][1]
    L1s = my[0]  if planedir=='X' else mx[0]  if planedir=='Y' else mx[0]
    L1e = my[-1] if planedir=='X' else mx[-1] if planedir=='Y' else mx[-1]
    L2s = mz[0]  if planedir=='X' else mz[0]  if planedir=='Y' else my[0]
    L2e = mz[-1] if planedir=='X' else mz[-1] if planedir=='Y' else my[-1]
    
    x = np.linspace(L1s,L1e,np.shape(V)[0])
    y = np.linspace(L2s,L2e,np.shape(V)[1])
    X,Y = np.meshgrid(x,y)
    X = X.T
    Y = Y.T
    
    #return V, X, Y, planedir, times
    return V, X, Y, planedir

###############################################################################

#---------- helper class for selectively hiding print statements
# https://stackoverflow.com/questions/8391411/suppress-calls-to-print-python

class _hide_print:
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._stdout

#-------------------------------------------------

def _get_varprops(lines: List[str]) -> List[Tuple[str, str, float]]:
    '''
    Pass in the list of lines represeting an FDS input file.
    Return list of tuples containing the variable name, direction ('X', 'Y', or 'Z'), 
        and location of the plane in the given direction.
    The first occurance of a variable name corresponds to plane 1.
    The second occurance of a variable name corresponds to plane 2, etc.
    (But note that a variable may be listed more than once for a given plane.)
    '''
    
    varprops=[]             # list of tuples(planeID, location), e.g. ('X', 12.5)
    varlines=''
    for line in lines:
        if line.find('&SLCF') == 0:
            
            #-------------- change e.g., QUANTITY='MASS FRACTION', SPEC_ID='OXYGEN' to
            #--------------     QUANTITY='OXYGEN MASS FRACTION'
            specid = re.search(r"SPEC_ID.*?=.*?'(.*?)'", line)
            if specid:
                specid = specid.group(1)
                quantity = re.search(r"QUANTITY.*?=.*?'(.*?)'",line)
                replacement = f"QUANTITY='{specid} {quantity.group(1)}'"
                line = re.sub(quantity.group(0), replacement, line)
                
            #-------------- append the line to all &SLCF lines
            
            varlines+=line
            
            if line.find('PBX') != -1:
                planeType = 'X'
                loc = float(re.search('PBX\s*=\s*(\d*\.?\d*)', line).group(1))
            elif line.find('PBY') != -1:
                planeType = 'Y'
                loc = float(re.search('PBY\s*=\s*(\d*\.?\d*)', line).group(1))
            elif line.find('PBZ') != -1:
                planeType = 'Z'
                loc = float(re.search('PBZ\s*=\s*(\d*\.?\d*)', line).group(1))
            else:
                sys.exit("ERROR: didn't find PBX, PBY, or PBZ on line: \n\t", line)
                
            varprops.append( (planeType, loc) )
                
            
            #-------------- add in vector components if present
            
            if re.search(r"VECTOR *=.TRUE.", line):
                varlines += "QUANTITY='U_VELOCITY', "
                varlines += "QUANTITY='V_VELOCITY', "
                varlines += "QUANTITY='W_VELOCITY', "
                varprops.append( (planeType, loc) )
                varprops.append( (planeType, loc) )
                varprops.append( (planeType, loc) )
                
    varlist = re.findall(r"QUANTITY.*?=.*?'(.*?)'", varlines)
    for i, var in enumerate(varlist):
        varprops[i] = (var, *varprops[i])
                
    return varprops

#-------------------------------------------------

def _get_mesh_topology(iflines: List[str]) -> Tuple[np.ndarray,np.ndarray,np.ndarray,int,int,int]:
    '''
    Pass in the list of lines in the input file.
    Return arrays: mx, my, mz are lists of mesh boundary locations in each direction,
        nmx, nmy, nmz are the number of meshes in each direction,
    Note the internal helper function get_nmxyz.
    '''
    
    #---------------------
    
    def _get_nmxyz(iflines: List[str], mult_id: str) -> Tuple[int, int, int]:
        '''
        Pass in the list of lines represeting an FDS input file.
        return tuple of numbers of meshes in x, y, z directions.
        '''
        ihi=ilo=jhi=jlo=khi=klo=0
        for line in iflines:
            if re.search('^&MULT', line) and re.search(f"ID\s*=\s*'{mult_id}'",line):
                
                i_lower = re.search('I_LOWER\s*=\s*(\d*)', line)
                if i_lower:
                    ilo = int(i_lower.group(1))
                i_upper = re.search('I_UPPER\s*=\s*(\d*)', line)
                if i_upper:
                    ihi = int(i_upper.group(1))
                    
                j_lower = re.search('J_LOWER\s*=\s*(\d*)', line)
                if j_lower:
                    jlo = int(j_lower.group(1))
                j_upper = re.search('J_UPPER\s*=\s*(\d*)', line)
                if j_upper:
                    jhi = int(j_upper.group(1))
                    
                k_lower = re.search('K_LOWER\s*=\s*(\d*)', line)
                if k_lower:
                    klo = int(k_lower.group(1))
                k_upper = re.search('K_UPPER\s*=\s*(\d*)', line)
                if k_upper:
                    khi = int(k_upper.group(1))
                    
        return ihi-ilo+1, jhi-jlo+1, khi-klo+1
    
    #---------------------
    
    for line in iflines:
        
        if re.search('^&MESH', line):
            xxyyzz = re.search('XB\s*=\s*(\d*\.?\d*\s*,\s*\d*\.?\d*,\s*\d*\.?\d*,\s*\d*\.?\d*,\s*\d*\.?\d*,\s*\d*\.?\d*)', line).group(1)
            xxyyzz = [float(i) for i in xxyyzz.split(',')]
            
            if not re.search('MULT_ID', line):
                mx = np.array([xxyyzz[0],xxyyzz[1]])
                mx = np.array([xxyyzz[2],xxyyzz[3]])
                mx = np.array([xxyyzz[4],xxyyzz[5]])
            else:
                mult_id = re.search("MULT_ID\s*=\s*'(.*)'", line).group(1)
                nmx,nmy,nmz = _get_nmxyz(iflines, mult_id)
                mx = xxyyzz[0] + (xxyyzz[1]-xxyyzz[0])*np.arange(0,nmx+1)
                my = xxyyzz[2] + (xxyyzz[3]-xxyyzz[2])*np.arange(0,nmy+1)
                mz = xxyyzz[4] + (xxyyzz[5]-xxyyzz[4])*np.arange(0,nmz+1)
            
            break
            
    return mx, my, mz, nmx, nmy, nmz
            

#-------------------------------------------------

def _get_meshes_ivarinmesh(varind: int, varprops: List[Tuple[str, str, float]],
                           mx: np.ndarray, my: np.ndarray, mz: np.ndarray) -> Tuple[List[int],List[int]]:
    '''
    Input varind: index of variable of interest in list of SLCF lines.
    Input varprops: list of tuples of variable properties e.g., (name, 'X', 15.0)
    Input mx, my, mz as the lists of mesh boundary locations in each direction.
    Return list of meshes that the plane intersects and index of the variable for that mesh.
    (If the plane is on a mesh boundary, the lower of the two meshes is considered in the returned arrays.)
    '''
    
    nmx = len(mx)-1      # number of meshes in x dir
    nmy = len(my)-1      # number of meshes in y dir
    nmz = len(mz)-1      # number of meshes in z dir
    
    mxyz = varprops[varind][1]
    xyz  = varprops[varind][2]
    
    meshes = []
    ivarinmesh = []
    if mxyz == 'Z':
        for k in range(nmz):
            if xyz >= mz[k] and xyz < mz[k+1]:
                break
        for j in range(nmy):
            for i in range(nmx):
                meshes.append( k*nmx*nmy + j*nmx + i)
                iii = -1
                for iv in range(varind+1):
                    if   varprops[iv][1]=='X' and varprops[iv][2] >= mx[i] and varprops[iv][2] <= mx[i+1]:
                        iii+=1
                    elif varprops[iv][1]=='Y' and varprops[iv][2] >= my[j] and varprops[iv][2] <= my[j+1]:
                        iii+=1
                    elif varprops[iv][1]=='Z' and varprops[iv][2] >= mz[k] and varprops[iv][2] <= mz[k+1]:
                        iii+=1
                ivarinmesh.append(iii)
                
    elif mxyz == 'Y':
        for j in range(nmy):
            if xyz >= my[j] and xyz < my[j+1]:
                break
        for k in range(nmz):
            for i in range(nmx):
                meshes.append( k*nmx*nmy + j*nmx + i)
                iii = -1
                for iv in range(varind+1):
                    if   varprops[iv][1]=='X' and varprops[iv][2] >= mx[i] and varprops[iv][2] <= mx[i+1]:
                        iii+=1
                    elif varprops[iv][1]=='Y' and varprops[iv][2] >= my[j] and varprops[iv][2] <= my[j+1]:
                        iii+=1
                    elif varprops[iv][1]=='Z' and varprops[iv][2] >= mz[k] and varprops[iv][2] <= mz[k+1]:
                        iii+=1
                ivarinmesh.append(iii)
                
    elif mxyz == 'X':
        for i in range(nmx):
            if xyz >= mx[i] and xyz < mx[i+1]:
                break
        for k in range(nmz):
            for j in range(nmy):
                meshes.append( k*nmx*nmy + j*nmx + i)
                iii = -1
                for iv in range(varind+1):
                    if   varprops[iv][1]=='X' and varprops[iv][2] >= mx[i] and varprops[iv][2] <= mx[i+1]:
                        iii+=1
                    elif varprops[iv][1]=='Y' and varprops[iv][2] >= my[j] and varprops[iv][2] <= my[j+1]:
                        iii+=1
                    elif varprops[iv][1]=='Z' and varprops[iv][2] >= mz[k] and varprops[iv][2] <= mz[k+1]:
                        iii+=1
                ivarinmesh.append(iii)
    else:
        sys.error("ERROR: value is", mxyz, "but must be 'X', or 'Y', or 'Z'")
        
    return meshes, ivarinmesh

#-------------------------------------------------

def _join_continued_lines(iflines: List[str]) -> List[str]:
    '''
    Input list of lines in the input file.
    Return list of relevant (non-empty or non-comment) lines
        multiline statements concatenated to a single line to aid processing.
    '''
    
    processed_lines = []
    nskip = 0
    for i,line in enumerate(iflines):
        if nskip > 0:
            nskip -= 1
            continue
        if re.search('^&', line):
            
            processed_lines.append(line)
            if not re.search('/', line):
                j=i
                while j<len(iflines):
                    j += 1
                    processed_lines[-1] = processed_lines[-1].rstrip() + ' ' + iflines[j]
                    if re.search('/', processed_lines[-1]):
                        nskip = j-i 
                        break
                        
    return processed_lines

###############################################################################

if __name__ == "__main__":
    
    inputFile = "case.fds"
    varName = "TEMPERATURE"
    pnumber = 2
    gskip = 1
    tskip = 1
    V, X, Y, pdir = get_var_slice(inputFile, varName, pnumber, gskip, tskip)
    
    xlabel = 'y (m)' if pdir=='X' else 'x (m)' if pdir=='Y' else 'x (m)'
    ylabel = 'z (m)' if pdir=='X' else 'z (m)' if pdir=='Y' else 'y (m)'
        
    #---------- plot results
    
    plt.rc('font', size=14)
    
    plt.figure(figsize=(8,6))
    plt.contourf(X,Y, V[:,:,-1], 80)
    plt.colorbar()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.axis('scaled')
    plt.title(varName);
    plt.show()
    
    nx,ny,nt = np.shape(V)
    Vavg = np.zeros((nx,ny))
    for it in range(nt):
        Vavg += V[:,:,it]
    Vavg /= nt
    
    plt.figure(figsize=(8,6))
    plt.contourf(X,Y,Vavg, 80)
    plt.colorbar()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.axis('scaled');    
    plt.title(f"Mean {varName}");
    plt.show()


