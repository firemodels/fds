import math as m
import numpy as np
import smop_util
# McGrattan
# 8-10-09
# reaction_rate.m

# input:
# T = temperature (K)
# Y = mass fraction

# output:
# r = reaction rate


def reaction_rate(T=None,Y=None):

    r=np.asarray([], dtype='object')
# temp/reaction_rate.m:14
    global dTdt
    global R0
    global E
    global A
    global residue
    r = smop_util.safe_set(r,(0,),np.dot(np.multiply(- A[0],Y[0]),m.exp(- E[0] / (np.multiply(R0,T)))) / dTdt)
# temp/reaction_rate.m:21
    r = smop_util.safe_set(r,(1,),np.dot(np.multiply(- A[1],Y[1]),m.exp(- E[1] / (np.multiply(R0,T)))) / dTdt)
# temp/reaction_rate.m:22
    r = smop_util.safe_set(r,(2,),np.dot(- residue[1],r[1]))
# temp/reaction_rate.m:23
    r=r.T
# temp/reaction_rate.m:24
    return r