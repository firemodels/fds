# This script calculates input parameters for an FDS simulation involving vegetative fuels.

# Modify the following lines for your particular type of fuel. 

Y_C = 0.497        # Mass fraction of carbon in the dry vegetation
Y_H = 0.059        # Mass fraction of hydrogen in the dry vegetation
Y_O = 0.406        # Mass fraction of oxygen in the dry vegetation
Y_C_char = 0.869   # Mass fraction of carbon in the char
Y_O_char = 0.049   # Mass fraction of oxygen in the char
nu_char = 0.16     # Fraction of mass of dry vegetation remaining as char after complete anaerobic pyrolysis
nu_ash = 0.0       # Fraction of mass of char remaining as ash after complete oxidation
E = 13100          # Energy released per unit mass oxygen consumed during combustion of pyrolyzate (kJ/kg)
W_pyr = 25         # Effective molecular weight of the pyrolyzate (g/mol)

# Calculated quantities. Some of these quantities are calculated automatically by FDS and are listed here as check.

Y_C_adj = Y_C*(1-nu_char*nu_ash)/(Y_C+Y_H+Y_O)
Y_H_adj = Y_H*(1-nu_char*nu_ash)/(Y_C+Y_H+Y_O)
Y_O_adj = Y_O*(1-nu_char*nu_ash)/(Y_C+Y_H+Y_O)

x = (Y_C_adj/12)/(Y_C_adj/12 + Y_H_adj + Y_O_adj/16 + nu_char*nu_ash)
y = (Y_H_adj   )/(Y_C_adj/12 + Y_H_adj + Y_O_adj/16 + nu_char*nu_ash)
z = (Y_O_adj/16)/(Y_C_adj/12 + Y_H_adj + Y_O_adj/16 + nu_char*nu_ash)

W_veg = (12*x+y+16*z)/(1-nu_char*nu_ash)
x_prime = nu_char*W_veg*(1-nu_ash)/(12*(1+Y_O_char/Y_C_char))
z_prime = 12*x_prime*Y_O_char/(16*Y_C_char)

nu_O2_char = (32*nu_char*W_veg*(1-nu_ash)-44*16*z_prime)/((44-32)*nu_char*W_veg)

nu_CO2 = (1+nu_O2_char-nu_ash)*(nu_char*W_veg)/44
nu_O2 = nu_O2_char*nu_char*W_veg/32
nu_pyr = (12*(x-x_prime)+y+16*(z-z_prime))/W_pyr

x_prime_prime = (x-x_prime)/nu_pyr
y_prime_prime = y/nu_pyr
z_prime_prime = (z-z_prime)/nu_pyr
W_pyr_check = 12*x_prime_prime+y_prime_prime+16*z_prime_prime

nu_O2_prime = (2*x_prime_prime+y_prime_prime/2-z_prime_prime)/2
Delta_h_pyr = (32*nu_O2_prime*E)/W_pyr

print()
print('The following values are used in the FDS input file:', end='\n\n')
print(f"{x_prime_prime:.3f}",' x_prime_prime, carbon subscript of pyrolyzate molecule (C on REAC or SPEC line)')
print(f"{y_prime_prime:.3f}",' y_prime_prime, hydrogen subscript of pyrolyzate molecule (H on REAC or SPEC line)')
print(f"{z_prime_prime:.3f}",' z_prime_prime, oxygen subscript of pyrolyzate molecule (O on REAC or SPEC line)')
print(f"{nu_O2_char:.3f}",' nu_O2_char, mass oxygen consumed per unit mass char oxidized (OXYGEN NU_SPEC on MATL line)')
print(f"{Delta_h_pyr:.0f}",' Delta_h_pyr, pyrolyzate heat of combustion (HEAT_OF_COMBUSTION on REAC line)')
print()
print('The following values are not used directly and can be used to check the calculation:', end='\n\n')
print(f"{Y_C_adj:.3f}",' Y_C_adj, adjusted carbon content')
print(f"{Y_H_adj:.3f}",' Y_H_adj, adjusted hydrogen content')
print(f"{Y_O_adj:.3f}",' Y_O_adj, adjusted oxygen content')
print(f"{x:.3f}",' x, carbon subscript for the dry vegetation molecule')
print(f"{y:.3f}",' y, hydrogen subscript for the dry vegetation molecule')
print(f"{z:.3f}",' z, oxygen subscript for the dry vegetation molecule')
print(f"{W_veg:.3f}",' W_veg, molecular weight of dry vegetation (g/mol)')
print(f"{x_prime:.3f}",' x_prime, carbon subscript in char molecule')
print(f"{z_prime:.3f}",' z_prime, oxygen subscript in char molecule')
print(f"{nu_CO2:.3f}",' nu_CO2, moles CO2 per moles char oxidized')
print(f"{nu_O2:.3f}",' nu_O2, moles O2 consumed per moles char oxidized')
print(f"{nu_pyr:.3f}",' nu_pyr, moles pyrolyzate per moles vegetation pyrolyzed')
print(f"{W_pyr_check:.2f}",' W_pyr_check, algebra check (g/mol)')
print(f"{nu_O2_prime:.3f}",' nu_O2_prime, moles O2 per mole pyrolyzate consumed')

