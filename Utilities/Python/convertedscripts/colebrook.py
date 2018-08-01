import numpy as np
# McDermott
# 2-23-09
# colebrook.m

# input:
# Re = Reynolds number
# RR = relative roughness
# ff = a guess for the value of f
# tol = convergence tolerance (relative error between f and ff)

# output:
# f = Moody friction factor based on Colebrook formula

# ref:
# Munson, Young and Okiishi, Fundamentals of Fluid Mechanics, Wiley, 1990.


def colebrook(Re=None,RR=None,ff=None,tol=None):

    iter_=0
# temp/colebrook.m:18
    error=tol + 1
# temp/colebrook.m:19
    while error > tol:

        iter_=iter_ + 1
# temp/colebrook.m:21
        f=(np.dot(- 2,np.log10(RR / 3.7 + 2.51 / (np.dot(Re,np.sqrt(ff)))))) ** (- 2)
# temp/colebrook.m:22
        error=np.fabs((f - ff) / ff)
# temp/colebrook.m:23
        ff=np.copy(f)
# temp/colebrook.m:24

    
    return f,error,iter_