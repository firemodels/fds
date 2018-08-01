import numpy as np
# McDermott
# 5-14-2009
# friction_factor_calc.m


def friction_factor_calc(dpdx=None,H=None,filename=None,mu_in=None):

    M=np.nan_to_num(np.genfromtxt(filename, delimiter = ',')[2:,0:], copy = False)
# temp/friction_factor_calc.m:7
    ubar=M[:,1]
# temp/friction_factor_calc.m:9
    if mu_in == None:
        mu=max(M[:,3])
# temp/friction_factor_calc.m:11
    else:
        mu=mu_in
# temp/friction_factor_calc.m:13
    
    rho=max(M[:,4])
# temp/friction_factor_calc.m:16
    U=ubar[ubar.size-1]
# temp/friction_factor_calc.m:18
    
    Re_H=np.dot(np.dot(H,U),rho) / mu
# temp/friction_factor_calc.m:19
    
    f_fds=np.dot(np.dot(2,(- dpdx)),H) / (np.dot(rho,U ** 2))
# temp/friction_factor_calc.m:20
    
    return f_fds,Re_H