
import numpy as np

R     = 8314.5
P     = 101325
W     = 28.8
mu    = 2.E-5
T_g_1 = 25 + 273.15
T_w_1 = 20 + 273.15
T_g_2 = 115 + 273.15
T_w_2 = 120 + 273.15
dx = 0.001
Vc = dx**3
L  = 0.01
Y_0 = 1.E-5 # initial mass fraction of aerosol
rhoa = 2000
ra = 1.E-6/2 # see input file
dTdx_1 = (T_g_1-T_w_1)/(0.5*dx)
dTdx_2 = (T_g_2-T_w_2)/(0.5*dx)
rho_g_1 = P*W/(R*T_g_1)
rho_g_2 = P*W/(R*T_g_2)

# FDS Tech Guide Eq. 8.9
alpha = 0.025/1.0 # see input file
Cs = 1.17
Ct = 2.2
Cm = 1.146
lam_1 = mu*np.sqrt(np.pi/(2*P*rho_g_1))
lam_2 = mu*np.sqrt(np.pi/(2*P*rho_g_2))
Kn_1 = lam_1/ra
Kn_2 = lam_2/ra
Cn_1 = 1 + 1.257*Kn_1 + 0.4*Kn_1*np.exp(-1.1/Kn_1)
Cn_2 = 1 + 1.257*Kn_2 + 0.4*Kn_2*np.exp(-1.1/Kn_2)
Kth_1 = 2*Cs*(alpha + Ct*Kn_1)*Cn_1/ ( (1 + 3*Cm*Kn_1)*(1 + 2*alpha + 2*Ct*Kn_1) )
Kth_2 = 2*Cs*(alpha + Ct*Kn_2)*Cn_2/ ( (1 + 3*Cm*Kn_2)*(1 + 2*alpha + 2*Ct*Kn_2) )

# note: for this case surface 1 is cold and soot deposits there
#       surface 2 is hot and soot moves away from it due to thermophoresis
u_th_1 = Kth_1*mu/rho_g_1/T_g_1 * dTdx_1; print(u_th_1)
u_th_2 = Kth_2*mu/rho_g_2/T_g_2 * dTdx_2; print(u_th_2)

T_END=2
f_1 = u_th_1*T_END / dx
f_2 = u_th_2*T_END / dx

# Y_0 = rhoa*Va/(rhoa*Va + rhog*(Vc-Va)); Y_0 = 1.E-5 per input file
Va_1 = Y_0*rho_g_1*Vc/( rhoa-Y_0*(rhoa-rho_g_1) )
Va_2 = Y_0*rho_g_2*Vc/( rhoa-Y_0*(rhoa-rho_g_2) )

ma_1 = rhoa*Va_1
ma_2 = rhoa*Va_2

Ma_1 = f_1*ma_1 *L**2/dx**2
Ma_2 = f_2*ma_2 *L**2/dx**2

# the value Ma_1 is the target mass deposited for this verification case
print(Ma_1,Ma_2)
