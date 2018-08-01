import numpy as np
import cmath as cm
# McDermott
# 4-8-10
# compression_wave_soln.m


def compression_wave_soln(rho0=None,x=None,y=None,a=None,c=None,t=None):

    b=np.sqrt(- 1 + a ** 2)
# temp/compression_wave_soln.m:7
    d=np.sqrt(- 1 + c ** 2)
# temp/compression_wave_soln.m:8
    x0=np.dot(2,np.arctan(np.dot(b / a,np.tan(np.arctan((1 + np.dot(a,np.tan(x / 2))) / b) - np.dot(b,t) / 2)) - 1 / a))
# temp/compression_wave_soln.m:10
    y0=np.dot(2,np.arctan(np.dot(d / c,np.tan(np.arctan((1 + np.dot(c,np.tan(y / 2))) / d) - np.dot(d,t) / 2)) - 1 / c))
# temp/compression_wave_soln.m:11
    Ix0=cm.log(- a ** 2 - np.cos(np.dot(2,np.arctan((1 + np.dot(a,np.tan(x0 / 2))) / b))) + np.dot(b,np.sin(np.dot(2,np.arctan((1 + np.dot(a,np.tan(x0 / 2))) / b)))))
# temp/compression_wave_soln.m:13
    Iy0=cm.log(- c ** 2 - np.cos(np.dot(2,np.arctan((1 + np.dot(c,np.tan(y0 / 2))) / d))) + np.dot(d,np.sin(np.dot(2,np.arctan((1 + np.dot(c,np.tan(y0 / 2))) / d)))))
# temp/compression_wave_soln.m:14
    Ix=cm.log(- a ** 2 - np.cos(np.dot(b,t) + np.dot(2,np.arctan((1 + np.dot(a,np.tan(x0 / 2))) / b))) + np.dot(b,np.sin(np.dot(b,t) + np.dot(2,np.arctan((1 + np.dot(a,np.tan(x0 / 2))) / b)))))
# temp/compression_wave_soln.m:16
    Iy=cm.log(- c ** 2 - np.cos(np.dot(d,t) + np.dot(2,np.arctan((1 + np.dot(c,np.tan(y0 / 2))) / d))) + np.dot(d,np.sin(np.dot(d,t) + np.dot(2,np.arctan((1 + np.dot(c,np.tan(y0 / 2))) / d)))))
# temp/compression_wave_soln.m:17
    q0=cm.log(rho0)
# temp/compression_wave_soln.m:19
    q=q0 + Ix - Ix0 + Iy - Iy0
# temp/compression_wave_soln.m:20
    rho=np.exp(q)
# temp/compression_wave_soln.m:22
    return rho