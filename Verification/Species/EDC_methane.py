
# Generate exact solution for EDC_methane verification cases

import numpy as np
from scipy.integrate import ode
import csv

def EDC_spec_methane(t, y):
    """
    Calculate derivatives for methane combustion system.

    Species indices:
    (0) - nitrogen
    (1) - methane
    (2) - oxygen
    (3) - co2
    (4) - h2o

    Args:
        t: Time variable
        y: Array of species concentrations

    Returns:
        dy: Array of derivatives for each species
    """

    m = 1                # chemical reaction coefficient fuel
    s = 2                # chemical reaction coefficient oxygen
    c = 1                # chemical reaction coefficient carbon dioxide
    h = 2                # chemical reaction coefficient water vapor
    tau_mix = 0.125      # set mixing time [s]

    # -----------------------
    # Logic for setting molar production rate
    # No suppression in this case
    # -----------------------
    r1 = min([y[1], y[2]/s])
    r_co2 = np.zeros(1)
    r_co2[0] = m * r1 / tau_mix

    dy = np.zeros(5)
    dy[0] = 0.0
    dy[1] = -(1/m) * r_co2[0]
    dy[2] = -(s/m) * r_co2[0]
    dy[3] = c * r_co2[0]
    dy[4] = h * r_co2[0]

    return dy


temp_0 = 298.15         # initial temperature [K]
volume = 0.001          # volume [m^3]
rho = 1.1294838         # density [kg/m^3]
R = 8.3145              # gas constant [J/mol K]
mass = rho * volume     # mass [kg]

yf0 = np.array([0.7247698, 0.0551642, 0.220066, 0.0, 0.0])  # initial mass fraction

y_MW = np.array([28.0134, 16.042460, 31.9988, 44.0095, 18.01528])  # [g/mol]

y_hf = np.array([0.0, -74873, 0.0, -393522, -241826])  # [J/mol]

N0 = (1000 * volume * rho * yf0) / y_MW  # initial moles

pres_0 = temp_0 * np.sum(N0) * R / volume  # initial pressure

tspan = np.arange(0, 2.1, 0.1)

y0 = np.array([N0[0], N0[1], N0[2], N0[3], N0[4]])

solver = ode(EDC_spec_methane)
solver.set_integrator('dopri5', rtol=1e-10, atol=1e-10)
solver.set_initial_value(y0, tspan[0])

# Integrate over time
T = [tspan[0]]
Y = [y0]
for t in tspan[1:]:
    solver.integrate(t)
    T.append(solver.t)
    Y.append(solver.y)

T = np.array(T)
Y = np.array(Y)

Nf = Y[-1, :]

# Convert back to mass fraction
yff = np.zeros((len(T), 5))
yff[:, 0] = (y_MW[0] * Y[:, 0]) / (1000 * rho * volume)
yff[:, 1] = (y_MW[1] * Y[:, 1]) / (1000 * rho * volume)
yff[:, 2] = (y_MW[2] * Y[:, 2]) / (1000 * rho * volume)
yff[:, 3] = (y_MW[3] * Y[:, 3]) / (1000 * rho * volume)
yff[:, 4] = (y_MW[4] * Y[:, 4]) / (1000 * rho * volume)

# Determine final temperature and pressure
tol0 = 1e5
tol = 0.1
Tf_guess = np.array([1000.0, 2000.0, 3000.0])  # final temperature guess [K]
count = 0

while abs(tol0) > tol:
    # ---------
    # NOTE: coefficients were found from NIST Webbook
    # ---------

    # initial temperature coeffs
    coeff0 = np.zeros((7, 5))
    coeff0[:, 0] = [28.98641, 1.853978, -9.647459, 16.63537, 0.000117, -8.671914, 0.0]
    coeff0[:, 1] = [-0.703029, 108.4773, -42.52157, 5.862788, 0.678565, -76.84376, -74.87310]
    coeff0[:, 2] = [31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 0.0]
    coeff0[:, 3] = [24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, -393.5224]
    coeff0[:, 4] = [30.09200, 6.832514, 6.793435, -2.534480, 0.082139, -250.8810, -241.8264]

    coeff = np.zeros((7, 5))

    # nitrogen cp coeffs [J/mol K]
    if Tf_guess[1] <= 500:
        coeff[:, 0] = [28.98641, 1.853978, -9.647459, 16.63537, 0.000117, -8.671914, 0.0]
    elif 500 < Tf_guess[1] <= 2000:
        coeff[:, 0] = [19.50583, 19.88705, -8.598535, 1.369784, 0.527601, -4.935202, 0.0]
    else:
        coeff[:, 0] = [35.51872, 1.128728, -0.196103, 0.014662, -4.553760, -18.97091, 0.0]

    # methane cp coeffs [J/mol K]
    if Tf_guess[1] <= 1300:
        coeff[:, 1] = [-0.703029, 108.4773, -42.52157, 5.862788, 0.678565, -76.84376, -74.87310]
    else:
        coeff[:, 1] = [85.81217, 11.26467, -2.114146, 0.138190, -26.42221, -153.5327, -74.87310]

    # oxygen cp coeffs [J/mol K]
    if Tf_guess[1] <= 700:
        coeff[:, 2] = [31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 0.0]
    elif 700 < Tf_guess[1] <= 2000:
        coeff[:, 2] = [30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 0.0]
    else:
        coeff[:, 2] = [20.91111, 10.72071, -2.020498, 0.146449, 9.245722, 5.337651, 0.0]

    # carbon dioxide cp coeffs [J/mol K]
    if Tf_guess[1] <= 1200:
        coeff[:, 3] = [24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, -393.5224]
    else:
        coeff[:, 3] = [58.16639, 2.720074, -0.492289, 0.038844, -6.447293, -425.9186, -393.5224]

    # water vapor cp coeffs [J/mol K]
    if Tf_guess[1] <= 1700:
        coeff[:, 4] = [30.09200, 6.832514, 6.793435, -2.534480, 0.082139, -250.8810, -241.8264]
    else:
        coeff[:, 4] = [41.96246, 8.622053, -1.499780, 0.098199, -11.15764, -272.1797, -241.8264]

    t = Tf_guess / 1000
    t0 = temp_0 / 1000
    del_h = np.zeros((5, 3))

    for j in range(3):
        for i in range(5):
            del_h[i, j] = ((coeff[0, i] * t[j] + (1/2) * coeff[1, i] * t[j]**2 + (1/3) * coeff[2, i] * t[j]**3 +
                           (1/4) * coeff[3, i] * t[j]**4 - coeff[4, i] / t[j] + coeff[5, i] - coeff[6, i]) -
                          (coeff0[0, i] * t0 + (1/2) * coeff0[1, i] * t0**2 + (1/3) * coeff0[2, i] * t0**3 +
                           (1/4) * coeff0[3, i] * t0**4 - coeff0[4, i] / t0 + coeff0[5, i] - coeff0[6, i])) * 1000

    h_fc = Nf[3] * y_hf[3] + Nf[4] * y_hf[4] - N0[1] * y_hf[1]
    del_hN = Nf[0] * del_h[0, :] + Nf[3] * del_h[3, :] + Nf[4] * del_h[4, :]
    RT0 = np.sum(N0) * R * temp_0
    RTf = -np.sum(Nf) * R * Tf_guess

    Tf_solve = h_fc + del_hN + RT0 + RTf
    tol0 = Tf_guess[0] - Tf_guess[2]

    if Tf_solve[1] < 0:
        if Tf_solve[0] > 0:
            mp = 0.5 * (Tf_guess[0] + Tf_guess[1])
            Tf_guess = np.array([Tf_guess[0], mp, Tf_guess[1]])
        else:
            mp = 0.5 * (Tf_guess[1] + Tf_guess[2])
            Tf_guess = np.array([Tf_guess[1], mp, Tf_guess[2]])
    else:
        if Tf_solve[0] < 0:
            mp = 0.5 * (Tf_guess[0] + Tf_guess[1])
            Tf_guess = np.array([Tf_guess[0], mp, Tf_guess[1]])
        else:
            mp = 0.5 * (Tf_guess[1] + Tf_guess[2])
            Tf_guess = np.array([Tf_guess[1], mp, Tf_guess[2]])

    count = count + 1

Tf = Tf_guess[1]
Pf = (Tf_guess[1] * np.sum(Nf) * R) / volume
dP = Pf - pres_0
Tf = Tf - 273.15

tss = np.array([1.0, 1.25, 1.5, 1.75, 2.0])
Tf_vec = np.array([Tf, Tf, Tf, Tf, Tf])
dP_vec = np.array([dP, dP, dP, dP, dP])

yff_output = np.column_stack((tspan, yff[:, 2], yff[:, 1], yff[:, 3], yff[:, 4]))

header1 = ['Time', 'CH4', 'O2', 'CO2', 'H2O']
with open('reactionrate_EDC_flim_1step_CH4_spec.csv', 'w', newline='') as fid:
    writer = csv.writer(fid)
    writer.writerow(header1)
    for j in range(len(tspan)):
        formatted_row = [f"{value:.6f}" for value in yff_output[j, :]]
        writer.writerow(formatted_row)

header1 = ['Time', 'TEMP', 'PRES']
with open('reactionrate_EDC_flim_1step_CH4_temppres.csv', 'w', newline='') as fid:
    writer = csv.writer(fid)
    writer.writerow(header1)
    for j in range(len(tss)):
        writer.writerow([f"{tss[j]:.6f}", f"{Tf_vec[j]:.6f}", f"{dP_vec[j]:.6f}"])




# Non premix case/fuel rich

temp_0 = 298.15        # initial temperature [K]
volume = 0.001         # volume [m^3]
rho = 1.1294838        # density [kg/m^3]
R = 8.3145             # gas constant [J/mol K]
mass = rho * volume    # mass [kg]

yf0 = np.array([0.7252, 0.0582, 0.2166, 0.0, 0.0])
y_MW = np.array([28.0134, 16.042460, 31.9988, 44.0095, 18.01528])  # [g/mol]

y_hf = np.array([0.0, -74873, 0.0, -393522, -241826])  # [J/mol]

N0 = (1000 * volume * rho * yf0) / y_MW  # initial moles

pres_0 = temp_0 * np.sum(N0) * R / volume  # initial pressure

tspan = np.arange(0, 65, 5)

y0 = np.array([N0[0], N0[1], N0[2], N0[3], N0[4]])

solver = ode(EDC_spec_methane)
#solver.set_integrator('dopri5', rtol=1e-10, atol=1e-10)
solver.set_integrator('vode', method='bdf')
solver.set_initial_value(y0, tspan[0])

Y = []
T = []
for t_val in tspan[1:]:
    if solver.successful():
        solver.integrate(t_val)
        Y.append(solver.y.copy())
        T.append(t_val)
    else:
        break

Y = np.array(Y)
T = np.array(T)

Nf = Y[-1, :]

# Convert back to mass fraction
yff = np.zeros((len(tspan), 5))
yff[1:, 0] = (y_MW[0] * Y[:, 0]) / (1000 * rho * volume)
yff[1:, 1] = (y_MW[1] * Y[:, 1]) / (1000 * rho * volume)
yff[1:, 2] = (y_MW[2] * Y[:, 2]) / (1000 * rho * volume)
yff[1:, 3] = (y_MW[3] * Y[:, 3]) / (1000 * rho * volume)
yff[1:, 4] = (y_MW[4] * Y[:, 4]) / (1000 * rho * volume)

# Determine final temperature and pressure
tol0 = 1e5
tol = 0.1
Tf_guess = np.array([1000.0, 2000.0, 3000.0])  # final temperature guess [K]
count = 0

while abs(tol0) > tol:
    # ---------
    # NOTE: coefficients were found from NIST Webbook
    # ---------

    # initial temperature coeffs
    coeff0 = np.zeros((7, 5))
    coeff0[:, 0] = np.array([28.98641, 1.853978, -9.647459, 16.63537, 0.000117, -8.671914, 0.0])
    coeff0[:, 1] = np.array([-0.703029, 108.4773, -42.52157, 5.862788, 0.678565, -76.84376, -74.87310])
    coeff0[:, 2] = np.array([31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 0.0])
    coeff0[:, 3] = np.array([24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, -393.5224])
    coeff0[:, 4] = np.array([30.09200, 6.832514, 6.793435, -2.534480, 0.082139, -250.8810, -241.8264])

    coeff = np.zeros((7, 5))

    # nitrogen cp coeffs [J/mol K]
    if Tf_guess[1] <= 500:
        coeff[:, 0] = np.array([28.98641, 1.853978, -9.647459, 16.63537, 0.000117, -8.671914, 0.0])
    elif 500 < Tf_guess[1] <= 2000:
        coeff[:, 0] = np.array([19.50583, 19.88705, -8.598535, 1.369784, 0.527601, -4.935202, 0.0])
    else:
        coeff[:, 0] = np.array([35.51872, 1.128728, -0.196103, 0.014662, -4.553760, -18.97091, 0.0])

    # methane cp coeffs [J/mol K]
    if Tf_guess[1] <= 1300:
        coeff[:, 1] = np.array([-0.703029, 108.4773, -42.52157, 5.862788, 0.678565, -76.84376, -74.87310])
    else:
        coeff[:, 1] = np.array([85.81217, 11.26467, -2.114146, 0.138190, -26.42221, -153.5327, -74.87310])

    # oxygen cp coeffs [J/mol K]
    if Tf_guess[1] <= 700:
        coeff[:, 2] = np.array([31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 0.0])
    elif 700 < Tf_guess[1] <= 2000:
        coeff[:, 2] = np.array([30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 0.0])
    else:
        coeff[:, 2] = np.array([20.91111, 10.72071, -2.020498, 0.146449, 9.245722, 5.337651, 0.0])

    # carbon dioxide cp coeffs [J/mol K]
    if Tf_guess[1] <= 1200:
        coeff[:, 3] = np.array([24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, -393.5224])
    else:
        coeff[:, 3] = np.array([58.16639, 2.720074, -0.492289, 0.038844, -6.447293, -425.9186, -393.5224])

    # water vapor cp coeffs [J/mol K]
    if Tf_guess[1] <= 1700:
        coeff[:, 4] = np.array([30.09200, 6.832514, 6.793435, -2.534480, 0.082139, -250.8810, -241.8264])
    else:
        coeff[:, 4] = np.array([41.96246, 8.622053, -1.499780, 0.098199, -11.15764, -272.1797, -241.8264])

    t = Tf_guess / 1000
    t0 = temp_0 / 1000

    del_h = np.zeros((5, 3))
    for j in range(3):
        for i in range(5):
            del_h[i, j] = ((coeff[0, i] * t[j] + (1/2) * coeff[1, i] * t[j]**2 + (1/3) * coeff[2, i] * t[j]**3 +
                           (1/4) * coeff[3, i] * t[j]**4 - coeff[4, i] / t[j] + coeff[5, i] - coeff[6, i]) -
                          (coeff0[0, i] * t0 + (1/2) * coeff0[1, i] * t0**2 + (1/3) * coeff0[2, i] * t0**3 +
                           (1/4) * coeff0[3, i] * t0**4 - coeff0[4, i] / t0 + coeff0[5, i] - coeff0[6, i])) * 1000

    h_fc = Nf[3] * y_hf[3] + Nf[4] * y_hf[4] - N0[1] * y_hf[1]
    del_hN = Nf[0] * del_h[0, :] + Nf[3] * del_h[3, :] + Nf[4] * del_h[4, :] + Nf[1] * del_h[1, :]
    RT0 = np.sum(N0) * R * temp_0
    RTf = -np.sum(Nf) * R * Tf_guess

    Tf_solve = h_fc + del_hN + RT0 + RTf
    tol0 = Tf_guess[0] - Tf_guess[2]

    if Tf_solve[1] < 0:
        if Tf_solve[0] > 0:
            mp = 0.5 * (Tf_guess[0] + Tf_guess[1])
            Tf_guess = np.array([Tf_guess[0], mp, Tf_guess[1]])
        else:
            mp = 0.5 * (Tf_guess[1] + Tf_guess[2])
            Tf_guess = np.array([Tf_guess[1], mp, Tf_guess[2]])
    else:
        if Tf_solve[0] < 0:
            mp = 0.5 * (Tf_guess[0] + Tf_guess[1])
            Tf_guess = np.array([Tf_guess[0], mp, Tf_guess[1]])
        else:
            mp = 0.5 * (Tf_guess[1] + Tf_guess[2])
            Tf_guess = np.array([Tf_guess[1], mp, Tf_guess[2]])

    count = count + 1

Tf = Tf_guess[1]
Pf = (Tf_guess[1] * np.sum(Nf) * R) / volume
dP = Pf - pres_0
Tf = Tf - 273.15

tss = 5 * (np.arange(1, 17)) + 25
Tf_vec = Tf * np.ones(16)
dP_vec = dP * np.ones(16)

tss = np.array([35, 40, 45, 50, 55, 60])
Tf_final = np.array([Tf, Tf, Tf, Tf, Tf, Tf])
dP_final = np.array([dP, dP, dP, dP, dP, dP])

yff[:, 0] = tspan

header1 = ['Time', 'CH4', 'O2', 'CO2', 'H2O']
with open('reactionrate_EDC_1step_CH4_nonmix_spec.csv', 'w', newline='') as fid:
    writer = csv.writer(fid)
    writer.writerow(header1)
    for j in range(7, len(tspan)):
        writer.writerow([f"{yff[j,0]:.6f}", f"{yff[j,1]:.6f}", f"{yff[j,2]:.6f}", f"{yff[j,3]:.6f}", f"{yff[j, 4]:.6f}"])

header1 = ['Time', 'TEMP', 'PRES']
with open('reactionrate_EDC_1step_CH4_nonmix_temppres.csv', 'w', newline='') as fid:
    writer = csv.writer(fid)
    writer.writerow(header1)
    for j in range(len(tss)):
        writer.writerow([f"{tss[j]:.6f}", f"{Tf_final[j]:.6f}", f"{dP_final[j]:.6f}"])


