"""

Script to produce exact numerical solution data for combined radiation and
conduction in slab between two blackbody surfaces at different temperatures.
(From M. Modest, Radiative Transfer, Second Edition, Sec. 21.2)

"""

import numpy as np
from scipy.special import expn
from scipy.sparse import diags
from matplotlib import pyplot as plt

def rhs( theta ):
    src = np.zeros(N_g - 2)
    for i in range( 0, len(src) ):
        for j in range(1, N_g):
            src[i] += (delta_tau/(2.*N))*(theta[j]**4 - theta[j-1]**4) \
                      *(expn(3, abs(i+1-j)*delta_tau) - expn(3, abs(i+2-j)*delta_tau))

    return src

# constants
sigma   = 5.67e-8           # W/m^2-K^4

# material properties
kappa_ar= [100., 2000.]     # m^-1
n       = 1.5
k       = 0.1               # W/m-K

# scenario parameters
T_1     = 700. + 273.15     # K
T_2     = 20. + 273.15      # K
L       = 0.1               # m
theta_L = T_2/T_1

# numerical paraters
N_g     = 500               # number of grid points
omega   = 0.1               # relaxation parameter
e2_c    = 1e-8              # convergence criterion

# initial guess for temperatures
x       = np.linspace(0., L, N_g)
T_con   = T_1 - (T_1 - T_2)*x/L
theta_init = T_con/T_1

# initialize solution array
T_e     = np.zeros((N_g, len(kappa_ar)))

# tridiagonal matrix for temperatures
A       = diags([1, -2, 1], [-1, 0, 1], shape=(N_g, N_g)).toarray()
A[0,0]  = 1.
A[0,1]  = 0.
A[-1,-1]= 1.
A[-1,-2]= 0.

# initialize right hand side nonlinear terms
b       = np.zeros( N_g )
b[0]    = 1.
b[-1]   = theta_L

# loop through absorption coeffcients
for i in range(0, len(kappa_ar)):

    # derived parameters
    N       = k*kappa_ar[i]/(4.*sigma*T_1**3)
    tau_L   = kappa_ar[i]*L
    delta_tau = tau_L/N_g

    # print derived parameters
    print("****************************")
    print("N            = ", N)
    print("tau_L        = ", tau_L)
    print("theta_L      = ", theta_L)
    print("delta_tau    = ", delta_tau)
    print("****************************")

    # solve nonlinear system of equations
    theta   = theta_init
    e2      = 1.
    while e2 > e2_c:
        b[1:N_g-1] = rhs( theta )
        theta_est = np.linalg.solve(A, b)
        e2 = np.sum( (theta_est - theta)**2 )/N_g
        print("Error = ", e2)
        theta = omega*theta_est + (1. - omega)*theta

    T_e[:,i] = theta*T_1 - 273.15

# format output data
output = np.hstack((x.reshape(-1,1), T_e))
head   = 'x (m), T (C) for kappa = ' + str(kappa_ar[0]) + ' (1/m), T (C) for kappa = '\
         + str(kappa_ar[1]) + ' (1/m)'

# save results to file
np.savetxt('radiation_conduction_exact.csv', output, delimiter=',', header=head)

