
# slab    Exact solution of the temperature inside a slab
#   with initial temperature T0 and both surfaces
#   subject to heat transfer from Tinf temperature from t = 0.
#
# T = slab(L,x,k,rho,cp,h,t,T0,Tinf,N)
#
#   L       Slab half-thickness (one side) (m)
#   x       distance from the slab center (m)
#   k       conductivity    (W/(K m))
#   rho     density         (kg/m3)
#   cp      heat capasity   (J/(kg K))
#   h       heat transfer coefficient (W/(m2 K))
#   t       time vector             (s)
#   T0      Initial temperature     (C)
#   Tinf    Gas temperature         (C)
#   N       Number of terms in the series
#   options (optional) options for fzero

import numpy as np
from scipy.optimize import brentq

def slab(L, x, k, rho, cp, h, t, T0, Tinf, N, options=None):
    if options is None:
        options = {}

    x = np.atleast_1d(x)
    t = np.atleast_1d(t)

    Bi = h * L / k
    alpha = k / (rho * cp)
    la = np.zeros(N)

    for n in range(1, N + 1):
        func = lambda val: (1.0 / np.tan(val * L)) - (val * L / Bi)
        lower_bound = (n - 1 + 0.00001) * np.pi / L
        upper_bound = (n - 0.00001) * np.pi / L
        xtol = options.get("TolX", 1e-12)
        la[n - 1] = brentq(func, lower_bound, upper_bound, xtol=xtol)

    laL = la * L
    Sla = np.zeros(N)
    for n in range(N):
        laN = laL[n]
        Sla[n] = np.sin(laN) / (laN + np.sin(laN) * np.cos(laN))

    theta0 = T0 - Tinf
    R = np.zeros((len(t), len(x)))

    for i in range(len(t)):
        if t[i] == 0:
            R[0, :] = np.ones(len(x))
        else:
            for n in range(N):
                S = Sla[n] * np.exp(-(la[n] ** 2) * alpha * t[i])
                for j in range(len(x)):
                    R[i, j] = R[i, j] + 2 * S * np.cos(la[n] * x[j])

    theta = R * theta0
    T = theta + Tinf
    return T

def main():

    # Parameters for Case A

    L = 0.1
    x = np.linspace(0, L, 6)
    k = 0.1
    rho = 100.
    cp = 1000.
    h = 100.
    t = np.array(range(0, 1801, 60), dtype=float)
    T0 = 20.
    Tinf = 120.
    N = 10

    temperatures = slab(L, x, k, rho, cp, h, t, T0, Tinf, N)

    output_file = "slab_temperatures.csv"

    # Stack time column with temperature matrix
    data_with_time = np.column_stack((t, temperatures))
    num_cols = data_with_time.shape[1]
    formats = ['%d'] + ['%.2f'] * (num_cols - 1)

    header = "Time,Back,8 cm,6 cm,4 cm,2 cm,Front"
    np.savetxt(output_file, data_with_time, delimiter=",",
               header=header, comments="", fmt=formats)
    print(f"Computed temperatures written to {output_file}")

if __name__ == "__main__":
    main()

