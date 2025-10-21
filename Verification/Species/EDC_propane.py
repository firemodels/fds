
# Verification of EDM for propane

import numpy as np
from scipy.integrate import ode
import csv

def EDC_spec_propane(t, y):
    """
    ODE system for propane combustion with EDC (Eddy Dissipation Concept)

    Species:
    (0) - nitrogen
    (1) - propane
    (2) - oxygen
    (3) - co2
    (4) - h2o
    """

    m = 1                # chemical reaction coefficient fuel
    s = 5                # chemical reaction coefficient oxygen
    c = 3                # chemical reaction coefficient carbon dioxide
    h = 4                # chemical reaction coefficient water vapor
    tau_mix = 0.125      # set mixing time [s]

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


temp_0 = 298.15        # initial temperature [K]
volume = 0.001         # volume [m^3]
rho = 1.2043424        # density [kg/m^3]
R = 8.3145             # gas constant [J/mol K]
mass = rho * volume    # mass [kg]

yf0 = np.array([0.720828, 0.060321, 0.218851, 0.0, 0.0])  # initial mass fraction
y_MW = np.array([28.0134, 44.09562, 31.9988, 44.0095, 18.01528])  # [g/mol]
y_hf = np.array([0.0, -103850, 0.0, -393522, -241826])  # [J/mol]

N0 = (1000 * volume * rho * yf0) / y_MW  # initial moles

pres_0 = temp_0 * np.sum(N0) * R / volume  # initial pressure

tspan = np.arange(0, 2.1, 0.1)
y0 = np.array([N0[0], N0[1], N0[2], N0[3], N0[4]])

# Setup ODE solver with tight tolerances
solver = ode(EDC_spec_propane)
solver.set_integrator('vode', method='adams', rtol=1e-10, atol=1e-10, max_step=0.1)
solver.set_initial_value(y0, tspan[0])

# Integrate over time span
Y = []
T = []
for t in tspan[1:]:
    if solver.successful():
        solver.integrate(t)
        Y.append(solver.y.copy())
        T.append(t)
    else:
        break

Y = np.array(Y)
T = np.array(T)

Nf = Y[-1, :]

# Convert back to mass fraction
yff = np.zeros((len(T), 6))
yff[:, 0] = T
yff[:, 1] = (y_MW[2] * Y[:, 2]) / (1000 * rho * volume)  # O2
yff[:, 2] = (y_MW[1] * Y[:, 1]) / (1000 * rho * volume)  # C3H8
yff[:, 3] = (y_MW[3] * Y[:, 3]) / (1000 * rho * volume)  # CO2
yff[:, 4] = (y_MW[4] * Y[:, 4]) / (1000 * rho * volume)  # H2O

# -----------------------
# Determine final temperature and pressure
# -----------------------
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
    coeff0[:, 1] = np.array([1, 1, 1, 1, 1, 1, 1])  # placeholder for propane
    coeff0[:, 2] = np.array([31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 0.0])
    coeff0[:, 3] = np.array([24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, -393.5224])
    coeff0[:, 4] = np.array([30.09200, 6.832514, 6.793435, -2.534480, 0.082139, -250.8810, -241.8264])

    # Initialize coefficient array
    coeff = np.zeros((7, 5))

    # nitrogen cp coeffs [J/mol K]
    if Tf_guess[1] <= 500:
        coeff[:, 0] = np.array([28.98641, 1.853978, -9.647459, 16.63537, 0.000117, -8.671914, 0.0])
    elif 500 < Tf_guess[1] <= 2000:
        coeff[:, 0] = np.array([19.50583, 19.88705, -8.598535, 1.369784, 0.527601, -4.935202, 0.0])
    else:
        coeff[:, 0] = np.array([35.51872, 1.128728, -0.196103, 0.014662, -4.553760, -18.97091, 0.0])

    # propane cp coeffs [J/mol K]
    if Tf_guess[1] <= 1300:
        coeff[:, 1] = np.array([1, 1, 1, 1, 1, 1, 1])  # placeholder for propane
    else:
        coeff[:, 1] = np.array([1, 1, 1, 1, 1, 1, 1])  # placeholder for propane

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
            del_h[i, j] = ((coeff[0, i] * t[j] + (1/2) * coeff[1, i] * t[j]**2 + (1/3) * coeff[2, i] * t[j]**3 + (1/4) * coeff[3, i] * t[j]**4 - coeff[4, i] / t[j] + coeff[5, i] - coeff[6, i]) \
                - (coeff0[0, i] * t0 + (1/2) * coeff0[1, i] * t0**2 + (1/3) * coeff0[2, i] * t0**3 + (1/4) * coeff0[3, i] * t0**4 - coeff0[4, i] / t0 + coeff0[5, i] - coeff0[6, i])) \
                * 1000

    # overwrite del_h[1, :] for propane (index 1 corresponds to propane)
    coeff2 = np.array([1.67536e-9, -5.46675e-6, 4.38029e-3, 2.7949, 6.11728e2])
    for j in range(3):
        del_h[1, j] = ((1/5) * coeff2[0] * Tf_guess[j]**5 + (1/4) * coeff2[1] * Tf_guess[j]**4 + (1/3) * coeff2[2] * Tf_guess[j]**3 + (1/2) * coeff2[3] * Tf_guess[j]**2 + coeff2[4] * Tf_guess[j]) \
            - ((1/5) * coeff2[0] * temp_0**5 + (1/4) * coeff2[1] * temp_0**4 + (1/3) * coeff2[2] * temp_0**3 + (1/2) * coeff2[3] * temp_0**2 + coeff2[4] * temp_0) \
            * (Tf_guess[j] - temp_0)

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

tss = np.array([1, 1.25, 1.5, 1.75, 2])
Tf_vec = np.array([Tf, Tf, Tf, Tf, Tf])
dP_vec = np.array([dP, dP, dP, dP, dP])

header1 = ['Time', 'O2', 'C3H8', 'CO2', 'H2O']
filename1 = '../../Verification/Species/reactionrate_EDC_flim_1step_C3H8_spec.csv'
with open('reactionrate_EDC_flim_1step_C3H8_spec.csv', 'w', newline='') as fid:
    writer = csv.writer(fid)
    writer.writerow(header1)
    for j in range(len(T)):
        writer.writerow([f"{yff[j,0]:.6f}", f"{yff[j,1]:.6f}", f"{yff[j,2]:.6f}", f"{yff[j,3]:.6f}", f"{yff[j,4]:.6f}"])

header1 = ['Time', 'TEMP', 'PRES']
with open('reactionrate_EDC_flim_1step_C3H8_temppres.csv', 'w', newline='') as fid:
    writer = csv.writer(fid)
    writer.writerow(header1)
    for j in range(len(tss)):
        writer.writerow([f"{tss[j]:.6f}", f"{Tf_vec[j]:.6f}", f"{dP_vec[j]:.6f}"])

