# C Weinschenk
# Verification of Mixture Fraction (Species,Temp,Pres)
# Fuel: Methane
# burke_schumann.m
# 9/2012
#
# Converted by Floyd
# 10-16-2025

import pandas as pd
import numpy as np
import os

# --- File Paths ---
base_dir = '../../Verification/Species/'
input_file_name = 'burke_schumann_devc.csv'

skip_case = False

if not os.path.exists(base_dir+input_file_name):
      skip_case = True
      print('Error: File ', input_file_name, ' does not exist. Skipping case.')

if skip_case: quit()

df = pd.read_csv(base_dir+input_file_name, skiprows=2, header=None)

burke = df.values

temperature = np.zeros((burke.shape[0], 15))
rho = np.zeros((burke.shape[0], 15))
h = np.zeros((burke.shape[0], 15))
hrrpuv = np.zeros((burke.shape[0], 15))
pressure = np.zeros((burke.shape[0], 15))
mix_frac = np.zeros((burke.shape[0], 15))
o2 = np.zeros((burke.shape[0], 15))
ch4 = np.zeros((burke.shape[0], 15))
h2o = np.zeros((burke.shape[0], 15))
co2 = np.zeros((burke.shape[0], 15))
n2 = np.zeros((burke.shape[0], 15))

for i in range(15):
   start_idx = 1 + 11 * i
   temperature[:, i] = burke[:, start_idx + 0] + 273.15
   rho[:, i] = burke[:, start_idx + 1]
   h[:, i] = burke[:, start_idx + 2]
   hrrpuv[:, i] = burke[:, start_idx + 3]
   pressure[:, i] = burke[:, start_idx + 4] + 101325
   mix_frac[:, i] = burke[:, start_idx + 5]
   o2[:, i] = burke[:, start_idx + 6]
   ch4[:, i] = burke[:, start_idx + 7]
   h2o[:, i] = burke[:, start_idx + 8]
   co2[:, i] = burke[:, start_idx + 9]
   n2[:, i] = burke[:, start_idx + 10]

n2_o2_ratio = n2[0, 0] / o2[0, 0] # Initial o2/n2 ratio

# Initialize as 2D arrays (time steps x 15 columns)
ox = np.zeros_like(o2)
prod = np.zeros_like(o2)

for i in range(o2.shape[1]):
   ox[:, i] = o2[:, i] + n2_o2_ratio * o2[:, i]
   prod[:, i] = h2o[:, i] + co2[:, i] + (n2[:, i] - n2_o2_ratio * o2[:, i])

# --- Calculate Expected Species and Temperature Profiles ---

volume = 0.001
temp_0 = temperature[0, 0]
density = rho[0, 4]
R = 8.314472

# Molecular weights and heats of formation
y_MW = np.array([28.0134, 16.042460, 31.9988, 44.0095, 18.015280]) # [g/mol]
y_hf = np.array([0.0, -74873, 0.0, 0.0, 0.0]) # [J/mol]

# Initial and final mass fractions (5 species x 15 probe points)
yf0 = np.vstack([n2[0, :], ch4[0, :], o2[0, :], co2[0, :], h2o[0, :]])
yff = np.vstack([n2[-1, :], ch4[-1, :], o2[-1, :], co2[-1, :], h2o[-1, :]]) # row end is index -1

# Initialize arrays for loop
N0 = np.zeros_like(yf0)
Nf = np.zeros_like(yff)
mean_mw_f = np.zeros(15)
Rw = np.zeros(15)

for i in range(15):
   y_MW_kg_mol = y_MW / 1000.0
   mass_point = volume * density # mass is constant for all 15 points in this block

   N0[:, i] = (mass_point * yf0[:, i]) / y_MW_kg_mol
   Nf[:, i] = (mass_point * yff[:, i]) / y_MW_kg_mol

   mole_frac_f = Nf[:, i] / np.sum(Nf[:, i])
   mean_mw_f[i] = np.sum(y_MW * mole_frac_f)

   Rw[i] = (R / mean_mw_f[i]) * 1000

# Heat of combustion [J/kg]
# hc = -y_hf(2)/y_MW(2)*1000; (y_hf(2) is index 1, y_MW(2) is index 1)
hc = -y_hf[1] / y_MW[1] * 1000
cp = 1000 # [J/kg/K] set constant for all species
cv = cp - Rw # [J/kg/K]

# --- Calculate State Relationships (Burke-Schumann) ---

z_st = yf0[1, 4]
f = mix_frac[0, :]

# Initialize arrays
Yf = np.zeros_like(f)
Yo2 = np.zeros_like(f)
Yp = np.zeros_like(f)
T_calc = np.zeros_like(f)

for i in range(len(f)):
   cv_i = cv[i]

   if f[i] > z_st and f[i] <= 1:
      # Fuel-rich side
      Yf[i] = (f[i] - z_st) / (1 - z_st)
      Yo2[i] = 0
      Yp[i] = (1 - f[i]) / (1 - z_st)
      T_calc[i] = f[i] * (-(z_st / (1 - z_st)) * (hc / cv_i)) + temp_0 + (z_st / ((1 - z_st) * cv_i)) * hc
   elif f[i] >= 0 and f[i] < z_st:
      # Fuel-lean side
      Yf[i] = 0
      Yo2[i] = 1 - f[i] / z_st
      Yp[i] = f[i] / z_st
      # T_calc(i) = f(i)*((hc)/(cv(i)))+temp_0;
      T_calc[i] = f[i] * (hc / cv_i) + temp_0
   else:
      # Stoichiometric point (should be covered by the previous cases, but for completeness)
      Yf[i] = 0
      Yo2[i] = 0
      Yp[i] = 1
      T_calc[i] = yf0[1, 4] * (hc / cv_i) + temp_0

T_adiab = T_calc[4]

# Normalized Temperatures
T_calc_norm = (T_calc - temp_0) / (T_adiab - temp_0)
# temperature(end,:) is the final time step's temperature row
FDS_temp_norm = (temperature[-1, :] - temp_0) / (T_adiab - temp_0)

# --- Write FDS and Expected Data to CSV Files ---

# 1. Write Expected Data
burke_expected = np.zeros((len(mix_frac[0, :]), 5)) # 15 rows, 5 columns

for i in range(len(mix_frac[0, :])):
    burke_expected[i, 0] = mix_frac[0, i]
    burke_expected[i, 1] = T_calc_norm[i]
    burke_expected[i, 2] = Yf[i]
    burke_expected[i, 3] = Yo2[i]
    burke_expected[i, 4] = Yp[i]

header1_expected = ['Mixture_Fraction', 'Temperature', 'Fuel', 'Ox', 'Prod']

# Create DataFrame for easy CSV writing with header
df_expected = pd.DataFrame(burke_expected, columns=header1_expected)
df_expected.to_csv(base_dir+'burke_schumann_expected.csv', index=False)

# 2. Write FDS Data
burke_FDS = np.zeros((len(mix_frac[0, :]), 5)) # 15 rows, 5 columns

for i in range(len(mix_frac[0, :])):
    burke_FDS[i, 0] = mix_frac[0, i]
    burke_FDS[i, 1] = FDS_temp_norm[i]
    burke_FDS[i, 2] = ch4[-1, i]
    burke_FDS[i, 3] = ox[-1, i]
    burke_FDS[i, 4] = prod[-1, i]

header1_FDS = ['Mixture_Fraction', 'FDS_Temperature', 'FDS_Fuel', 'FDS_Ox', 'FDS_Prod']

# Create DataFrame for easy CSV writing with header
df_FDS = pd.DataFrame(burke_FDS, columns=header1_FDS)
df_FDS.to_csv(base_dir+'burke_schumann_FDS.csv', index=False)
