"""
McDermott
7-21-11
wavelet_transform.py

Python version of WAVELET_ERROR function in dump.f90
"""

import numpy as np
import matplotlib.pyplot as plt

# signal
S = np.array([0, 1, 0, 1], dtype=float)

# plot signal
plt.plot(np.arange(1, len(S) + 1), S, 'bo--')
plt.xlabel('cell index')
plt.ylabel('normalized signal')
plt.show()

# normalize signal
SMAX = np.max(S)
SMIN = np.min(S)
DS = SMAX - SMIN

if DS < 1e-6:
    WAVELET_ERROR = 0.0
else:
    SS = (S - SMIN) / DS

    # discrete Haar wavelet transform
    M = 2
    N = M
    A = np.zeros((M, N))
    C = np.zeros((M, N))

    for I in range(M):
        for J in range(N):
            K = 2 * J
            if I == 0:
                A[I, J] = 0.5 * (SS[K] + SS[K + 1])
                C[I, J] = 0.5 * (SS[K] - SS[K + 1])
            else:
                A[I, J] = 0.5 * (A[I - 1, K] + A[I - 1, K + 1])
                C[I, J] = 0.5 * (A[I - 1, K] - A[I - 1, K + 1])
        N //= 2

    C1 = np.sum(C[0, :])
    C2 = np.sum(C[1, :])

    WAVELET_ERROR = abs(C1 - C2)

    print(f"WAVELET_ERROR = {WAVELET_ERROR:.6f}")
    print("A =\n", A)
    print("C =\n", C)


