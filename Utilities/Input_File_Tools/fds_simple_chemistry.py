# Floyd
# fds_simple_chemsitry.m
# 10/22/2025

# Replicate the basic procesing of a simple chemistry REAC line ouptutting
# lumped species definitions, reaction stoichiometry, EPUMO2, and
# constant gamma specific heat.

import math
import numpy as np

#### The following section is input data from the REAC line
# Define the fuel C_m H_n O_a N_b
m = 3
n = 8
a = 0
b = 0
# Set HoC to -1 to use EPUMO2 and vice versa
HoC = -1 # kJ/kg
Ideal = False
EPUMO2 = 13100 # kJ/kg

# post flame CO, Soot, and HCN yields (kg/kg of fuel)
y_CO = 0.012
y_Soot = 0.024
Soot_H_fraction = 0 # fraction of atoms in Soot that are H
y_HCN = 0.0

# 1 or 2 step
chem_steps = 1

# 2-step parameters
fuel_C_to_CO = 2./3.
fuel_H_to_H2 = 0.
fuel_N_to_HCN = 0.

# Soot tracked separately
soot_separate = False

# Ambient air
y_O2_Infty = 0.232378
y_CO2_Infty = 0.000595
T_amb = 20
P_amb = 101325
P_STP = 101325
Humidity = 40 # %

# Water Evap Prop
T_ref_H2O = 273.15
h_v_ref_H2O = 2500792.3981106
T_boil = 373.15

#Other constants
R0 = 8314.472
gamma = 1.4

if (Soot_H_fraction < 0 or Soot_H_fraction > 1):
   print('Soot_H_fraction must be >= 0 and <= 1')
   quit()

if (chem_steps != 1 and chem_steps !=2):
   print('chem_steps must be 1 or 2')
   quit()

if (m < 0 or n < 0 or a < 0 or b < 0):
   print('Fuel m,n,a,b values must be >=0')
   quit()

if (HoC > 0 and EPUMO2 > 0):
   print('Only one of HoC or EPUMO2 can be > 0')
   quit()

atom_weight ={'H':1.00794,'C':12.0107,'N':14.0067,'O':15.9994}
spec_MW = {
   'Fuel':m*atom_weight['C']+n*atom_weight['H']+a*atom_weight['O']+b*atom_weight['N'],
   'N2':2*atom_weight['N'],
   'O2':2*atom_weight['O'],
   'CO2':atom_weight['C']+2*atom_weight['O'],
   'CO':atom_weight['C']+atom_weight['O'],
   'Soot':(1-Soot_H_fraction)*atom_weight['C']+Soot_H_fraction*atom_weight['H'],
   'H2':2*atom_weight['H'],
   'H2O':2*atom_weight['H']+atom_weight['O'],
   'HCN':atom_weight['H']+atom_weight['C']+atom_weight['N']
   }

if (a==0 and b==0):
   print('Fuel is C',m,'H',n,'with MW=',f'{spec_MW['Fuel']:.5f}','g/mol')
elif (a==0 and b!=0):
   print('Fuel is C',m,'H',n,'N',b,'with MW=',spec_MW['Fuel'],'g/mol')
else:
   print('Fuel is C',m,'H',n,'O',a,'N',b,'with MW=',spec_MW['Fuel'],'g/mol')

spec_H_F = {
   'Fuel':-1,
   'N2':0,
   'O2':0,
   'CO2':-393.51,
   'CO':-110.535196,
   'Soot':0,
   'H2':0,
   'H2O':-241.826,
   'HCN':133.08246
   } #kJ/mol

dry_air_MW = 1/((1-y_O2_Infty-y_CO2_Infty)/spec_MW['N2']+y_O2_Infty/spec_MW['O2']+y_CO2_Infty/spec_MW['CO2'])

def cp_g_H2O_f(T):
   cp_g_H2O = -0.394796083E+05/T**2+0.575573102E+03/T+0.931782653E+00+0.722271300E-02*T+\
   -0.734256000E-05*T**2+0.495504000E-08*T**3-0.133693000E-11*T**4
   cp_g_H2O = cp_g_H2O*R0/spec_MW['H2O']
   return cp_g_H2O

def cp_f_H2O_f(T):
# Compute water vapor gas enthalpy
   cp_f_H2O = 1326371304./T**2-24482953.88/T+187942.8776-767.899505*T+\
   1.761556813*T**2-0.002151167128*T**3+1.092570813E-06*T**4
   cp_f_H2O = cp_f_H2O*R0/spec_MW['H2O']
   return cp_f_H2O

h_g_H2O = np.zeros(375)
h_f_H2O = np.zeros(375)

for i in range(len(h_f_H2O)):
   cp_g_2 = cp_g_H2O_f(max(200,i))
   cp_f_2 = cp_f_H2O_f(max(200,i))
   if (i>0):
      h_g_H2O[i] = h_g_H2O[i-1]+0.5*(cp_g_1+cp_g_2)
      h_f_H2O[i] = h_f_H2O[i-1]+0.5*(cp_f_1+cp_f_2)
   cp_g_1 = cp_g_2
   cp_f_1 = cp_f_2

h_f = spec_H_F['H2O']*1000/spec_MW['H2O']*1000
h_g_ref = h_g_H2O[298] + 0.15*(h_g_H2O[299]-h_g_H2O[298])
h_g_H2O = h_g_H2O + (h_f - h_g_ref )
i_T = int(T_ref_H2O)
h_g_ref = h_g_H2O[i_T] + (T_ref_H2O-i_T)*(h_g_H2O[i_T+1]-h_g_H2O[i_T])
h_f_ref = h_f_H2O[i_T] + (T_ref_H2O-i_T)*(h_f_H2O[i_T+1]-h_f_H2O[i_T])
h_v = h_g_ref - h_f_ref
h_f_H2O = h_f_H2O - (h_v_ref_H2O - h_v)
h_v_H2O = h_g_H2O - h_f_H2O

def x_H2O_f(T,Humidity):
   # returns ambient water vapor mole fraction
   pratio = P_amb/P_STP
   i_T = int(T)
   h_v =  h_v_H2O[i_T] + (T-i_T)*(h_v_H2O[i_T+1]-h_v_H2O[i_T])
   DHOR_T = h_v*spec_MW['H2O']/R0
   T_Boil_Eff = DHOR_T*T_boil/(DHOR_T-T_boil*math.log(pratio))
   i_T = int(T_Boil_Eff)
   h_v_TBE = h_v_H2O[i_T] + (T_Boil_Eff-i_T)*(h_v_H2O[i_T+1]-h_v_H2O[i_T])
   DHOR_TBE = h_v_TBE*spec_MW['H2O']/R0
   DHOR = 0.5*(DHOR_T+DHOR_TBE)
   x_H2O =  min(1.,math.exp(DHOR*(1/T_Boil_Eff-1/T)))*Humidity/100
   return x_H2O

def y_H2O_f(TC,Humidity):
   # returns ambient water vapor mass fraction
   T = TC + 273.15
   x_H2O = x_H2O_f(T,Humidity)
   y_H2O = x_H2O/(dry_air_MW/spec_MW['H2O']+(1-dry_air_MW/spec_MW['H2O'])*x_H2O)
   return y_H2O

species = ['Fuel','N2','O2','CO2','CO','Soot','H2','H2O','HCN']

i_Fuel = 0
i_N2 = 1
i_O2 = 2
i_CO2 = 3
i_CO = 4
i_Soot = 5
i_H2 = 6
i_H2O = 7
i_HCN = 8

y_H2O_amb = y_H2O_f(T_amb,Humidity)

# Spec list is per spec_MW
y_air = np.array([0,(1-y_O2_Infty-y_CO2_Infty)*(1-y_H2O_amb),y_O2_Infty*(1-y_H2O_amb),y_CO2_Infty*(1-y_H2O_amb),0,0,0,y_H2O_amb,0])
x_air = np.zeros(9)
i = 0
for key in spec_MW:
   x_air[i] = y_air[i]/spec_MW[key]
   i+=1
x_air = x_air/sum(x_air)

i = 0
h_f_air = 0
for key in spec_MW:
   h_f_air = h_f_air + x_air[i] * spec_H_F[key]
   i+=1

i=0
mw_air = 0
for key in spec_MW:
   mw_air = mw_air + x_air[i]*spec_MW[key]
   i+=1

print('')
print('Ambient conditions are:')
print('Temperature:',f'{T_amb:.2f}','C')
print('Pressure:',f'{P_amb:.0f}','Pa')
print('Humidity:',f'{Humidity:.2f}','%')

y_fuel = np.array([1,0,0,0,0,0,0,0,0])
x_fuel = y_fuel

nu_reac_T = np.zeros(9)
yield_T = np.zeros(9)
y_prod = np.zeros(9)
x_prod = np.zeros(9)

nu_reac_T[i_Fuel]=-1
nu_reac_T[i_CO] = y_CO * spec_MW['Fuel']/spec_MW['CO']
nu_reac_T[i_Soot] = y_Soot * spec_MW['Fuel']/spec_MW['Soot']
nu_reac_T[i_HCN] = y_HCN * spec_MW['Fuel']/spec_MW['HCN']
nu_reac_T[i_CO2] = m - nu_reac_T[i_CO] - nu_reac_T[i_Soot] * (1-Soot_H_fraction) - nu_reac_T[i_HCN]
if (nu_reac_T[i_CO2] < 0):
   print('nu_CO2 is ',f'{nu_reac_T[i_CO2]:.4f}','. Check fuel defintion and post-flame yields.')
nu_reac_T[i_H2O] = (n -  nu_reac_T[i_Soot] * Soot_H_fraction - nu_reac_T[i_HCN])/2
if (nu_reac_T[i_H2O] < 0):
   print('nu_H2O is ',f'{nu_reac_T[i_H2O]:.4f}','. Check fuel defintion and post-flame yields.')
nu_reac_T[i_N2] = (b - nu_reac_T[i_HCN])/2
if (nu_reac_T[i_N2] < 0):
   print('nu_N2 is ',f'{nu_reac_T[i_N2]:.4f}','. Check fuel defintion and post-flame yields.')
nu_reac_T[i_O2] = -(nu_reac_T[i_CO2]+(nu_reac_T[i_CO]+nu_reac_T[i_H2O]-a)/2)

x_prod[i_N2] = -x_air[i_N2] / x_air[i_O2] * nu_reac_T[i_O2] + nu_reac_T[i_N2]
x_prod[i_CO2] = - x_air[i_CO2] / x_air[i_O2] * nu_reac_T[i_O2] + nu_reac_T[i_CO2]
x_prod[i_H2O] = -x_air[i_H2O]  / x_air[i_O2] * nu_reac_T[i_O2] + nu_reac_T[i_H2O]
x_prod[i_CO] = nu_reac_T[i_CO]
x_prod[i_Soot] = nu_reac_T[i_Soot]
x_prod[i_HCN] = nu_reac_T[i_HCN]

if (soot_separate):
   x_prod = x_prod /(sum(x_prod) - x_prod[i_Soot])
else:
   x_prod = x_prod / sum(x_prod)

i = 0
h_f_prod = 0
for key in spec_MW:
   h_f_prod = h_f_prod + x_prod[i] * spec_H_F[key]
   i+=1

i = 0
for key in spec_MW:
   if (not soot_separate or (soot_separate and key!='Soot')):
      y_prod[i] = x_prod[i]*spec_MW[key]
   i+=1
i=0
mw_prod = 0
for key in spec_MW:
   if (not soot_separate or (soot_separate and key!='Soot')):
      mw_prod = mw_prod + x_prod[i]*spec_MW[key]
   i+=1

y_prod = y_prod/sum(y_prod)

print('')
if (HoC > 0):
   if (Ideal):
      iHoC = HoC
      rHoC = (HoC*spec_MW['Fuel']/1000+nu_reac_T[i_CO]*(spec_H_F['CO2']-spec_H_F['CO'])+ \
                                       nu_reac_T[i_Soot]*(spec_H_F['CO2']-spec_H_F['Soot'])+ \
                                       nu_reac_T[i_HCN]*(spec_H_F['CO2']-spec_H_F['HCN'])+ \
                                       nu_reac_T[i_HCN]*(spec_H_F['H2O']/2-spec_H_F['HCN']))/spec_MW['Fuel']*1000
   else:
      rHoC = HoC
      iHoC = (rHoC*spec_MW['Fuel']/1000-nu_reac_T[i_CO]*(spec_H_F['CO2']-spec_H_F['CO'])- \
                                        nu_reac_T[i_Soot]*(spec_H_F['CO2']-spec_H_F['Soot'])- \
                                        nu_reac_T[i_HCN]*(spec_H_F['CO2']-spec_H_F['HCN'])- \
                                        nu_reac_T[i_HCN]*(spec_H_F['H2O']/2-spec_H_F['HCN']))/spec_MW['Fuel']*1000

   EPUMO2 = rHoC*nu_reac_T[i_Fuel]/nu_reac_T[i_O2]*spec_MW['Fuel']/spec_MW['O2']

else:
   rHoC = EPUMO2*nu_reac_T[i_O2]/nu_reac_T[i_Fuel]*spec_MW['O2']/spec_MW['Fuel']
   iHoC = (rHoC*spec_MW['Fuel']/1000-nu_reac_T[i_CO]*(spec_H_F['CO2']-spec_H_F['CO'])- \
                                     nu_reac_T[i_Soot]*(spec_H_F['CO2']-spec_H_F['Soot'])- \
                                     nu_reac_T[i_HCN]*(spec_H_F['CO2']-spec_H_F['HCN'])- \
                                     nu_reac_T[i_HCN]*(spec_H_F['H2O']/2-spec_H_F['HCN']))/spec_MW['Fuel']*1000

spec_H_F['Fuel'] = iHoC*spec_MW['Fuel']/1000+m*spec_H_F['CO2']+n/2*spec_H_F['H2O']
nu_fuel = -1
m_fuel = -1
nu_air = nu_reac_T[i_O2]/x_air[i_O2]
m_air = -m_fuel/spec_MW['Fuel']*nu_air*mw_air
m_prod = -(m_fuel+m_air)
if (soot_separate):
   m_prod_s = y_Soot
   m_prod = m_prod - y_Soot
   nu_prod_s = nu_reac_T[i_Soot]
   nu_prod = -(nu_fuel*spec_MW['Fuel']+nu_air*mw_air -nu_reac_T[i_Soot]*spec_MW['Soot'] )/mw_prod
else:
   nu_prod = -(nu_fuel*spec_MW['Fuel']+nu_air*mw_air)/mw_prod

fuel_C_to_CO = 2./3.
fuel_H_to_H2 = 0.
fuel_N_to_HCN = 0.

if (chem_steps==2):
   nu_reac_1 = np.zeros(9)
   nu_reac_2 = np.zeros(9)
   y_int = np.zeros(9)
   x_int = np.zeros(9)
   nu_reac_1[i_Fuel]=-1
   nu_reac_1[i_CO] = fuel_C_to_CO * m
   nu_reac_1[i_H2] = fuel_H_to_H2 * n/2
   nu_reac_1[i_HCN] = fuel_N_to_HCN * b
   nu_reac_1[i_Soot] = (m - nu_reac_1[i_CO] - nu_reac_1[i_HCN])/(1-Soot_H_fraction)
   nu_reac_1[i_N2] = (1-fuel_N_to_HCN) * b/2
   nu_reac_1[i_H2O] = n/2 - nu_reac_1[i_H2] - (nu_reac_1[i_Soot]*Soot_H_fraction+nu_reac_1[i_HCN])/2
   nu_reac_1[i_O2] = -(nu_reac_1[i_CO]+nu_reac_1[i_H2O]-a)/2
   x_int[i_N2] = -x_air[i_N2] / x_air[i_O2] * nu_reac_1[i_O2] + nu_reac_1[i_N2]
   x_int[i_H2O] = -x_air[i_H2O]  / x_air[i_O2] * nu_reac_1[i_O2] + nu_reac_1[i_H2O]
   x_int[i_CO2]= -x_air[i_CO2]  / x_air[i_O2] * nu_reac_1[i_O2]
   x_int[i_CO] = nu_reac_1[i_CO]
   x_int[i_Soot] = nu_reac_1[i_Soot]
   x_int[i_HCN] = nu_reac_1[i_HCN]
   s_x_int = sum(x_int)
   s_x_int_2 = s_x_int
   if (soot_separate):
      s_x_int = s_x_int - x_int[i_Soot]
   x_int = x_int / s_x_int

   i = 0
   h_f_int = 0
   for key in spec_MW:
      h_f_int = h_f_int + x_int[i] * spec_H_F[key]
      i+=1

   HoC1 = (spec_H_F['Fuel'] - nu_reac_1[i_CO]*spec_H_F['CO'] -\
                              nu_reac_1[i_HCN]*spec_H_F['HCN'] -\
                              nu_reac_1[i_H2O]*spec_H_F['H2O'])/spec_MW['Fuel']*1000
   EPUMO2_1 = HoC1*nu_reac_1[i_Fuel]/nu_reac_1[i_O2]*spec_MW['Fuel']/spec_MW['O2']
   fuel_H_F_1 = iHoC*spec_MW['Fuel']/1000+nu_reac_1[i_CO]*spec_H_F['CO']

   nu_reac_2 = nu_reac_T - nu_reac_1
   nu_reac_2[i_Fuel] = -1
   nu_reac_2[i_O2] = -(nu_reac_2[i_CO2]+(nu_reac_2[i_CO]+nu_reac_2[i_H2O])/2)
   nu_reac_2 = nu_reac_2 / s_x_int_2

   i = 0
   for key in spec_MW:
      if (not soot_separate or (soot_separate and key!='Soot')):
         y_int[i] = x_int[i]*spec_MW[key]
      i+=1
   y_int = y_int/sum(y_int)

   i=0
   mw_int = 0
   for key in spec_MW:
      if (not soot_separate or (soot_separate and key!='Soot')):
         mw_int = mw_int + x_int[i]*spec_MW[key]
      i+=1

   nu_air_1 = nu_reac_1[i_O2]/x_air[i_O2]
   m_air_1 = -m_fuel/spec_MW['Fuel']*nu_air_1*mw_air
   m_int_1 = -(m_fuel+m_air_1)
   if (soot_separate):
      nu_int_1_s = nu_reac_1[i_Soot]
      m_int_1_s = nu_reac_1[i_Soot]*spec_MW['Soot']/spec_MW['Fuel']
      m_int_1 = m_int_1 - m_int_1_s
      nu_int_1 = -(nu_fuel*spec_MW['Fuel']+nu_air_1*mw_air+nu_reac_1[i_Soot]*spec_MW['Soot'])/mw_int
   else:
      nu_int_1 = -(nu_fuel*spec_MW['Fuel']+nu_air_1*mw_air)/mw_int

   nu_air_2 = nu_reac_2 [i_O2]/x_air[i_O2]
   m_air_2 = -m_fuel/mw_int*nu_air_2*mw_air
   m_int_2 = -(m_fuel+m_air_2)
   nu_int_2 = -(nu_fuel*mw_int+nu_air_2*mw_air)/mw_prod
   if (soot_separate):
      nu_int_2_s = nu_reac_2[i_Soot]
      m_int_2_s = nu_reac_2[i_Soot]*spec_MW['Soot']/spec_MW['Fuel']
      m_int_2 = m_int_2 - m_int_2_s
      nu_int_2 = -(nu_fuel*mw_int+nu_air_2*mw_air+nu_int_2_s*spec_MW['Soot'])/mw_prod
   else:
      nu_int_2 = -(nu_fuel*mw_int+nu_air_2*mw_air)/mw_prod

   HoC2 = (h_f_int - h_f_air * nu_air_2 - nu_int_2*h_f_prod)/mw_int*1000
   EPUMO2_2 = -HoC2/nu_reac_2[i_O2]*mw_int/spec_MW['O2']

   print('')
   print('Intermediate Reactions Stoichiometry')
   print('Species\t\t  Nu R1  \t  Nu R2')
   for i in range(len(species)):
      print(species[i],'\t\t',f'{nu_reac_1[i]:.6f}','\t',f'{nu_reac_2[i]:.6f}')

   print('')
   print('Step 1 HoC:        ',f'{HoC1:.2f}',' kJ/kg')
   print('Step 1 EPUMO2      ',f'{EPUMO2_1:.2f}',' kJ/kg')
   print('')
   print('Reaction 1 Lumped Species Stoich Coeff.')
   print('Species\t\t  Moles \t  Mass')
   print('Fuel\t\t',f'{nu_fuel:.5f}','\t',f'{m_fuel:.5f}')
   print('Air\t\t',f'{nu_air_1:.5f}','\t',f'{m_air_1:.5f}')
   print('Int. Prod.\t',f'{nu_int_1:.5f}','\t',f'{m_int_1:.5f}')
   if (soot_separate):
         print('Soot\t\t',f'{nu_int_1_s:.5f}','\t',f'{m_int_1_s:.5f}')
   print('')
   print('Step 2 HoC:        ',f'{HoC2:.2f}',' kJ/kg')
   print('Step 2 EPUMO2      ',f'{EPUMO2_2:.2f}',' kJ/kg')
   print('')
   print('Reaction 2 Lumped Species Stoich Coeff.')
   print('Species\t\t  Moles \t  Mass')
   print('Int. Prod.\t',f'{nu_fuel:.5f}','\t',f'{m_fuel:.5f}')
   print('Air\t\t',f'{nu_air_2:.5f}','\t',f'{m_air_2:.5f}')
   print('Products\t',f'{nu_int_2:.5f}','\t',f'{m_int_2:.5f}')
   if (soot_separate):
         print('Soot\t\t',f'{nu_int_2_s:.5f}','\t',f'{m_int_2_s:.5f}')

print('')
print('Total Reaction Stoichiometry')
print('Species\t\t   Nu   \t  Yield')
for i in range(len(species)):
   yield_T[i] = -nu_reac_T[i]*spec_MW[species[i]]/(nu_reac_T[i_Fuel]*spec_MW['Fuel'])
   print(species[i],'\t\t',f'{nu_reac_T[i]:.6f}','\t',f'{yield_T[i]:.6f}')

print('')
print('Total Reaction Lumped Species Stoich Coeff.')
print('Species\t\t  Moles \t  Mass')
print('Fuel\t\t',f'{nu_fuel:.5f}','\t',f'{m_fuel:.5f}')
print('Air\t\t',f'{nu_air:.5f}','\t',f'{m_air:.5f}')
print('Prodcuts\t',f'{nu_prod:.5f}','\t',f'{m_prod:.5f}')
if (soot_separate):
   print('Soot\t\t',f'{nu_prod_s:.5f}','\t',f'{m_prod_s:.5f}')

print('')
print('Fuel H_F:    ',f'{spec_H_F['Fuel']:.3f}',' kJ/mol')
print('Ideal HoC:   ',f'{iHoC:.2f}',' kJ/kg')
print('Reazlied HoC:',f'{rHoC:.2f}',' kJ/kg')
print('Fuel EPUMO2: ',f'{EPUMO2:.2f}',' kJ/kg')
if (EPUMO2 < 11000 or EPUMO2 > 16000):
   print('EPUMO2 is not close to 13100 kJ/kg, check fuel molecule definition and HoC setting')

if (chem_steps==1):
   print('')
   print('Mole Fractions for Lumped Species')
   print('Species\t   Air  \t  Fuel \t\tProducts')
   for i in range(len(species)):
      if (not soot_separate or (soot_separate and i!=i_Soot)):
         print(species[i],'\t',f'{x_air[i]:.6f}','\t',f'{x_fuel[i]:.6f}','\t',f'{x_prod[i]:.6f}')

   print('')
   print('Mass Fractions for Lumped Species')
   print('Species\t   Air  \t  Fuel \t\tProducts')
   for i in range(len(species)):
      if (not soot_separate or (soot_separate and i!=i_Soot)):
         print(species[i],'\t',f'{y_air[i]:.6f}','\t',f'{y_fuel[i]:.6f}','\t',f'{y_prod[i]:.6f}')

   print('')
   print('Lumped Species Properties')
   print('Species\t\t MW(g/mol)\tHeat of Formation (kJ/mol)')
   print('Air\t\t',f'{mw_air:.2f}','\t\t',f'{h_f_air:.4f}')
   print('Fuel\t\t',f'{spec_MW['Fuel']:.2f}','\t\t',f'{spec_H_F['Fuel']:.4f}')
   print('Products\t',f'{mw_prod:.2f}','\t\t',f'{h_f_prod:.4f}')
   if (soot_separate):
      print('Soot\t\t',f'{spec_MW['Soot']:.2f}','\t\t',f'{spec_H_F['Soot']:.4f}')

else: #chem_steps==2
   print('')
   print('Mole Fractions for Lumped Species')
   print('Species\t   Air  \t  Fuel \t\tInt Products\tProducts')
   for i in range(len(species)):
      if (not soot_separate or (soot_separate and i!=i_Soot)):
         print(species[i],'\t',f'{x_air[i]:.6f}','\t',f'{x_fuel[i]:.6f}','\t',f'{x_int[i]:.6f}','\t',f'{x_prod[i]:.6f}')

   print('')
   print('Mass Fractions for Lumped Species')
   print('Species\t   Air  \t  Fuel \t\tInt Products\tProducts')
   for i in range(len(species)):
      if (not soot_separate or (soot_separate and i!=i_Soot)):
         print(species[i],'\t',f'{y_air[i]:.6f}','\t',f'{y_fuel[i]:.6f}','\t',f'{y_int[i]:.6f}','\t',f'{y_prod[i]:.6f}')

   print('')
   print('Lumped Species Properties')
   print('Species\t\t MW(g/mol)\tHeat of Formation (kJ/mol)')
   print('Air\t\t',f'{mw_air:.4f}','\t',f'{h_f_air:.4f}')
   print('Fuel\t\t',f'{spec_MW['Fuel']:.4f}','\t',f'{spec_H_F['Fuel']:.4f}')
   print('Int. Prod.\t',f'{mw_int:.4f}','\t',f'{h_f_int:.4f}')
   print('Products\t',f'{mw_prod:.4f}','\t',f'{h_f_prod:.4f}')
   if (soot_separate):
      print('Soot\t\t',f'{spec_MW['Soot']:.4f}','\t',f'{spec_H_F['Soot']:.4f}')

