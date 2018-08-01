import math as m
import numpy as np
# Antonellis
# 06-23-2010
# section2_soln.m


def section2_soln(rho0=None,x=None,y=None,B=None,w=None,t=None):

    x0=np.dot(2,m.atan(np.dot(m.tan(x / 2),m.exp(np.dot(- B / w,m.sin(np.dot(w,t)))))))
# temp/section2_soln.m:7
    y0=np.dot(2,m.atan(np.dot(m.tan(y / 2),m.exp(np.dot(- B / w,m.sin(np.dot(w,t)))))))
# temp/section2_soln.m:8
    q0=np.log(rho0)
# temp/section2_soln.m:10
    q=q0 + np.log((1 + np.dot((m.tan(x0 / 2)) ** 2.0,m.exp(np.dot(np.dot(2,B) / w,m.sin(np.dot(w,t)))))) / (1 + (m.tan(x0 / 2)) ** 2)) + np.log((1 + np.dot((m.tan(y0 / 2)) ** 2.0,m.exp(np.dot(np.dot(2,B) / w,m.sin(np.dot(w,t)))))) / (1 + (m.tan(y0 / 2)) ** 2)) - np.dot(np.dot(2,B) / w,m.sin(np.dot(w,t)))
# temp/section2_soln.m:11
    
    
    rho=m.exp(q)
# temp/section2_soln.m:15
    return rho