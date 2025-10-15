# McDermott
# 11-1-2016
# water_evap_1_const_gamma.py
#
# Convert to python by Floyd
# 10/15/2025
#
# Calculation of expected values for water_evap_1_const_gamma.fds verification case.
# See fds/Verification/Sprinklers_and_Sprays/water_evap_1_const_gamma.csv

import os

# --- Constants ---
H_L = -1.9757398E+004 # J/kg (Latent Heat of Vaporization from problem setup)
R = 8314.472 # Pa * m^3 / (kmol * K) (Universal Gas Constant)
gam = 1.4 # Specific heat ratio (gamma)
P_1 = 101325 # Pa (Initial Pressure)
V = 1 # m^3 (Volume)
T_1 = 200 + 273.15 # K (Initial Gas Temperature)
T_w = 20 + 273.15 # K (Water/Evaporation Temperature - Note: T_w is not used in the final calculation, but is included for completeness)
M_w = 0.01 # kg (Mass of Water evaporated)

W_a = 28.84834 # kg/kmol (Molecular Weight of Air)
W_w = 18.01528 # kg/kmol (Molecular Weight of Water)

# --- Ideal Gas Specific Heats (J/(kg*K)) ---
# cv = R/W * 1/(gam-1)
cv_a = R/W_a * 1 / (gam - 1)
cv_w = R/W_w * 1 / (gam - 1)

# cp = R/W * gam/(gam-1)
cp_a = R/W_a * gam / (gam - 1)
cp_w = R/W_w * gam / (gam - 1)

# --- Compute Initial State ---
# Ideal gas law: RHO = P * W / (R * T)
RHO_1 = P_1 * W_a / (R * T_1)
# Mass = RHO * V
M_a = RHO_1 * V # Initial mass of air

# --- Final State Calculation ---

# Determine final mass fraction of water vapor in gas
Y_w = M_w / (M_a + M_w)
Y_a = 1 - Y_w

# Compute new mixture molecular weight
# W_mix = 1 / (Y_w/W_w + Y_a/W_a)
W = 1 / (Y_w / W_w + Y_a / W_a)

# Compute final density
RHO_2 = (M_w + M_a) / V

# Determine final temperature (T_2) from energy balance (Internal energy is conserved)
# M_a*cv_a*T_1 = M_a*cv_a*T_2 + M_w*cv_w*T_2 - H_L (where H_L is negative)
# T_2 = ( M_a*cv_a*T_1 + H_L ) / (M_a*cv_a + M_w*cv_w)
T_2 = (M_a * cv_a * T_1 + H_L) / (M_a * cv_a + M_w * cv_w)

# Determine final pressure (P_2) from ideal gas law for the mixture
# P = RHO * R * T / W_mix
P_2 = RHO_2 * R * T_2 / W

# --- Changes in Properties ---

# Change in pressure, Pa
dP = P_2 - P_1

# Change in gas enthalpy, J (dH = M_mix * cp_mix * T_2 - M_a * cp_a * T_1)
# Since $c_{p,mix} = Y_a c_{p,a} + Y_w c_{p,w}$ and $M_{mix} = M_a + M_w$, this is:
# $dH = (M_a c_{p,a} + M_w c_{p,w}) T_2 - M_a c_{p,a} T_1$
# NOTE: The formula in the original MATLAB code seems to assume that the water vapor
# is not present in the initial state and contributes to the final enthalpy.
# This assumes the water enters the gas phase at $T_2$.
dH = (M_a * cp_a + M_w * cp_w) * T_2 - M_a * cp_a * T_1

# --- Output Results to CSV ---

# Define output directory path
# This assumes the script is run from a location where the relative path works
# In a real-world scenario, you might want a more robust path handling (e.g., using pathlib)
ddir = '../../Verification/Sprinklers_and_Sprays/'
output_filename = ddir+'water_evap_1_const_gamma.csv'

# Ensure the directory exists (optional, but good practice)
os.makedirs(os.path.dirname(output_filename), exist_ok=True)

# Open and write to the file
with open(output_filename, 'w') as fid:
    # Header row
    header = ['Time', 'Rel. Hum', 'h_gas', 'h_water', 'dens', 'temp', 'pres', 'vapor']
    fid.write(','.join(header) + '\n')

    # Data for time 8 and 10
    # Note: MATLAB's fprintf formatting is used for precision matching
    # dH and H_L are converted to kJ (divided by 1000)
    # T_2 is converted to Celsius (T_2 - 273.15)
    # RHO_2 is delta RHO (RHO_2 - RHO_1)
    # M_w is M_w

    # Format strings match: %i, %11.9f, %14.7e, %14.7e, %6.4f, %10.7f, %14.7e, %6.4f
    data_format = "{:d},{:11.9f},{:14.7e},{:14.7e},{:6.4f},{:10.7f},{:14.7e},{:6.4f}\n"

    # Line 1: Time = 8
    line_8 = data_format.format(
        8,
        2.109675497,
        dH / 1000,
        -H_L / 1000,
        RHO_2 - RHO_1,
        T_2 - 273.15,
        dP,
        M_w
    )
    fid.write(line_8)

    # Line 2: Time = 10
    line_10 = data_format.format(
        10,
        2.109675497,
        dH / 1000,
        -H_L / 1000,
        RHO_2 - RHO_1,
        T_2 - 273.15,
        dP,
        M_w
    )
    fid.write(line_10)

print(f"Calculation complete. Results written to: {output_filename}")