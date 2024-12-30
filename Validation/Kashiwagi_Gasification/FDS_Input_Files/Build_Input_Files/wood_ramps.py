"""
McDermott
7-18-2024
wood_ramps.py

Generate thermal property ramps from Lautenberger and Pello CNF 2009 paper, Table 1.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

T_r   = 300.  # K
T_max = 5300. # K
T   = np.array([20., 300., 600., 900., 1200.]) + 273.
SIGMA=5.670373E-8 # W/m2/K4

# dry wood

c_0 = 1.664 # kJ/kg K
n_c = 0.660
c = c_0*(T/T_r)**n_c

k_0 = 0.176 # W/m K
n_k = 0.594
k = k_0*(T/T_r)**n_k

print('')
print('dry wood properties:')

cdf = pd.DataFrame({'T': T-273., 'c': c})
# print(cdf)
T_values = T
F_values = c
print('')
for TT, FF in zip(T_values, F_values):
    print(f"&RAMP ID='c_v dry wood', T={TT}, F={FF:.4g}/")

kdf = pd.DataFrame({'T': T-273., 'k': k})
# print(kdf)
T_values = T
F_values = k
print('')
for TT, FF in zip(T_values, F_values):
    print(f"&RAMP ID='k dry wood', T={TT}, F={FF:.4g}/")

# plt.figure(1)
# plt.plot(T,c)

# plt.figure(2)
# plt.plot(T,k)

# plt.show()

# char

c_0 = 1.219 # kJ/kg K
n_c = 0.283
c = c_0*(T/T_r)**n_c

k_0 = 0.065 # W/m K
n_k = 0.435
gamma = 3.3e-3
k = k_0*(T/T_r)**n_k + gamma*SIGMA*T**3

print('')
print('char properties:')

cdf = pd.DataFrame({'T': T-273., 'c': c})
# print(cdf)
T_values = T
F_values = c
print('')
for TT, FF in zip(T_values, F_values):
    print(f"&RAMP ID='c_v char', T={TT}, F={FF:.4g}/")

kdf = pd.DataFrame({'T': T-273., 'k': k})
# print(kdf)
T_values = T
F_values = k
print('')
for TT, FF in zip(T_values, F_values):
    print(f"&RAMP ID='k + k_r char', T={TT}, F={FF:.4g}/")

# ash

c_0 = 1.244 # kJ/kg K
n_c = 0.315
c = c_0*(T/T_r)**n_c

k_0 = 0.058 # W/m K
n_k = 0.353
gamma = 6.4e-3
k = k_0*(T/T_r)**n_k + gamma*SIGMA*T**3

print('')
print('ash properties:')

cdf = pd.DataFrame({'T': T-273., 'c': c})
# print(cdf)
T_values = T
F_values = c
print('')
for TT, FF in zip(T_values, F_values):
    print(f"&RAMP ID='c_v ash', T={TT}, F={FF:.4g}/")

kdf = pd.DataFrame({'T': T-273., 'k': k})
# print(kdf)
T_values = T
F_values = k
print('')
for TT, FF in zip(T_values, F_values):
    print(f"&RAMP ID='k + k_r ash', T={TT}, F={FF:.4g}/")




