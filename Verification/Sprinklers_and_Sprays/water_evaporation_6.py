m_w   = 1
c_w   = 1
m_g_0 = 1
c_p_g = 1
m_d   = 0.01
R     = 8.314
W     = 28
m_g_f = m_g_0 + m_d
T_g_0 = 373.15
c_l   = 2
T_b   = 373.15
T_d_0 = 293.15
h_v   = 173.15
T_w_0 = 523.15
V     = 1

T_g_f = (-m_g_0*c_p_g*T_g_0+m_d*(T_b-T_d_0)-c_p_g*T_b*m_d+m_d*h_v+R*m_g_0*T_g_0/W-m_w*c_w*T_w_0) / (-m_w*c_w-m_g_0*c_p_g-m_d*c_p_g+R*m_g_f/W)
Delta_p = R*(m_g_f*T_g_f - m_g_0*T_g_0)/(W*V)
print('The final temperature is ' + str(T_g_f-273.15) + ' C')
print('The pressure rise is ' + str(Delta_p) + ' kPa')
