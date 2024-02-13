L       = 1.76
b_f     = 0.22
t_f     = 0.005
t_w     = 0.006
t_i     = 0.01
d       = 0.58
rho_a   = 2700.
rho_i   = 20.
V_a     = L*(b_f*t_f*2+(d-2*t_f)*t_w)
V_c     = 0.8
V_i     = (L+2*t_i)*( (b_f+2*t_i)*(t_f+2*t_i)*2 + (d+2*t_i-2*(t_f+2*t_i))*(t_w+2*t_i) ) - V_a
p_0     = 101.325
R       = 0.008314
W       = 0.028
T_g_0   = 1073.15
rho_g   = p_0*W/(R*T_g_0)
c_g     = 1.
c_a     = 0.9
c_i     = 1.
T_a_0   = 293.15
V_g_1   = V_c - V_a
n_1     = p_0*V_g_1/(R*T_g_0)
U_g     = rho_g*V_g_1*c_g*T_g_0 - p_0*V_g_1
U_a     = rho_a*V_a*c_a*T_a_0
U_i     = rho_i*V_i*c_i*T_a_0
T_1     = (U_g+U_a)/(rho_a*V_a*c_a + rho_g*V_g_1*c_g - n_1*R)
V_g_2   = V_c - V_a - V_i
n_2     = p_0*V_g_2/(R*T_g_0)
T_2     = (U_g+U_a+U_i)/(rho_i*V_i*c_i + rho_a*V_a*c_a + rho_g*V_g_2*c_g - n_2*R)
print('n_1 is ' + str(n_1) + ' mol')
print('n_2 is ' + str(n_2) + ' mol')
print('U_g is ' + str(U_g) + ' kJ')
print('U_a is ' + str(U_a) + ' kJ')
print('U_i is ' + str(U_i) + ' kJ')
print('T_1 is ' + str(T_1-273.15) + ' C')
print('T_2 is ' + str(T_2-273.15) + ' C')

