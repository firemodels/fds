
# This script just computes the exact solution of the flat_fire case.

import numpy as np

V_0 = 400
h = 8
g = 9.8
C_d = 0.2
rho_a = 1.2
rho_d = 1000
D = 5e-3

K = 3 * rho_a * C_d / (4 * rho_d * D)

dt = 0.05
tend = 1.65

tvec = np.arange(0, tend + dt, dt)

xexact = np.log(V_0 * K * tvec + 1) / K
yexact = h \
    + (g / (2 * (K * V_0)**2)) * np.log(V_0 * K * tvec + 1) \
    - g * tvec ** 2 / 4 \
    - g * tvec / (2 * K * V_0)
uexact = V_0 / (V_0 * K * tvec + 1)
vexact = g / (2 * K * V_0 * (K * V_0 * tvec + 1)) \
    - g * tvec / 2 \
    - g / (2 * K * V_0)

mat = np.column_stack((tvec, xexact, yexact, uexact, vexact))
np.savetxt('flat_fire.csv', mat, delimiter=',', fmt=['%.2f', '%.5f', '%.5f', '%.5f', '%.5f'])

