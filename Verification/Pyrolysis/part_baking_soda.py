
# Generate exact solutions for part_baking_soda cases

import numpy as np
import pandas as pd

rho = 2200
R = 8.3145
A = 3.4e11
E = 103000
r_0 = 2.5e-6
T = np.array([420,450,500])
k = A * np.exp(-E / (R * T))
pd.options.display.float_format = '{:.4f}'.format

for i in range(len(T)):
    t = np.linspace(0, 10, 101)
    r1 = (r_0**3 * np.exp(-k[i]*t))**(1/3)  # first-order model
    rs = r_0 * (1 - k[i] * t)  # spherical contraction model
    d1 = np.maximum(0,r1 * 2e6)
    ds = np.maximum(0,rs * 2e6)
    Ts = str(int(T[i]))

    df = pd.DataFrame({ 'Time (s)': t, f'Dia (mu-m) First-order {Ts} K': d1, f'Dia (mu-m) Spherical {Ts} K': ds })
    df.to_csv(f'part_baking_soda_{Ts}K.csv', index=False, float_format='%.4f')

