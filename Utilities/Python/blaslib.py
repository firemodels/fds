#McDermott
#03-05-2016
#blaslib.py

import math
import numpy as np

def blasius_analytic(u0,zmax,mu,rho,xc):
    # print 'test'
    # return

    # etamax = maximum eta for calculation
    # steps = number of steps between 0 and etamax
    # fppwall = initial (wall) value of 2nd derivative of Blasius function
    # outputs:
    # eta - the similarity coordinate normal to the wall
    # f, fp, fpp, fppp - the Blasius function and it first 3 derivatives

    etamax=zmax/math.sqrt(mu/rho*xc/u0)
    steps=257
    deta = etamax/(steps-1)
    eta = np.zeros(steps)
    f = np.zeros(steps)
    fp = np.zeros(steps)
    fpp = np.zeros(steps)
    fppp = np.zeros(steps)
    k1 = np.zeros(3)
    k2 = np.zeros(3)
    k3 = np.zeros(3)
    k4 = np.zeros(3)

    # initial guess for fpp
    fpp[0] = 0.3318

    for i in range(steps-1):

        eta[i+1] = eta[i] + deta

        # predictor
        #1st
        k1[0]=fp[i]
        k1[1]=fpp[i]
        k1[2]=-f[i]*fpp[i]/2.

        fbar = f[i] + 0.5*deta * k1[0]
        fpbar = fp[i] + 0.5*deta * k1[1]
        fppbar = fpp[i] + 0.5*deta * k1[2]

        #2nd
        k2[0]=fpbar
        k2[1]=fppbar
        k2[2]=-fbar*fppbar/2.

        fbar = f[i] + 0.5*deta * k2[0]
        fpbar = fp[i] + 0.5*deta * k2[1]
        fppbar = fpp[i] + 0.5*deta * k2[2]

        #3rd
        k3[0]=fpbar
        k3[1]=fppbar
        k3[2]=-fbar*fppbar/2.
        
        fbar = f[i] + deta * k3[0]
        fpbar = fp[i] + deta * k3[1]
        fppbar = fpp[i] + deta * k3[2]

        #4th
        k4[0]=fpbar
        k4[1]=fppbar
        k4[2]=-fbar*fppbar/2.

        # corrector
        f[i+1] = f[i] + deta * (k1[0]+2.*k2[0]+2.*k3[0]+k4[0])/6.
        fp[i+1] = fp[i] + deta * (k1[1]+2.*k2[1]+2.*k3[1]+k4[1])/6.
        fpp[i+1] = fpp[i] + deta * (k1[2]+2.*k2[2]+2.*k3[2]+k4[2])/6.
        fppp[i+1] = -f[i+1]*fpp[i+1]/2.

    return eta, fp











