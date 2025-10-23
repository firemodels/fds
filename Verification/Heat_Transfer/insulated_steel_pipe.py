
# Write out exact solution of the heat conduction equation in radial coordinates

import math

h = 0.0
k = [0.0] * 4 
r = [0.0] * 5 
T = [0.0] * 6 
q = 0.0
pi = 0.0
TT = 0.0
rr = 0.0
i = 0
iz = 0
pi = math.atan(1.0) * 4.0

# Set parameters
h = 10.0
k[1] = 0.2
k[2] = 50.0
k[3] = 0.2
r[1] = 0.01
r[2] = r[1] + 0.02
r[3] = r[2] + 0.01
r[4] = r[3] + 0.02
T[0] = 20.0
T[5] = 480.0

# Calculate heat flux q
q = (T[5] - T[0]) / (
    (1.0 / (2 * pi * r[1] * h)) +
    math.log(r[2] / r[1]) / (2 * pi * k[1]) +
    math.log(r[3] / r[2]) / (2 * pi * k[2]) +
    math.log(r[4] / r[3]) / (2 * pi * k[3]) +
    (1.0 / (2 * pi * r[4] * h))
)

# Calculate interface temperatures
T[1] = T[0] + q / (2 * pi * r[1] * h)
T[2] = T[1] + math.log(r[2] / r[1]) * q / (2 * pi * k[1])
T[3] = T[2] + math.log(r[3] / r[2]) * q / (2 * pi * k[2])
T[4] = T[3] + math.log(r[4] / r[3]) * q / (2 * pi * k[3])

# Open file for writing (unit 10 in Fortran)
with open('insulated_steel_pipe.csv', 'w') as f10:
    # Write header
    f10.write('Radius,Temp\n')

    # Loop from 100 to 0 (inclusive) with step -1
    for i in range(100, -1, -1):
        rr = r[1] + i * (r[4] - r[1]) / 100.0

        # Determine layer index iz
        if rr < r[2]:
            iz = 2
        if r[2] <= rr and rr < r[3]:
            iz = 3
        if r[3] <= rr:
            iz = 4

        # Calculate temperature at depth
        TT = T[iz] + math.log(rr / r[iz]) * (T[iz - 1] - T[iz]) / math.log(r[iz - 1] / r[iz])

        # Write data line (format: f7.4, comma, f7.2)
        f10.write(f'{rr:7.4f},{TT:7.2f}\n')


