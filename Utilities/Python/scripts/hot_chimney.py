
def hot_chimney():
    import csv
    from pathlib import Path

    import numpy as np
    
    # Given parameters
    W = 0.18 # m
    D = 0.84 # m
    H = 4.0  # m
    Tw = 30.0 + 273.15 # K (Wall temperature)
    T_inf = 20.0 + 273.15 # K (Ambient temperature)
    g = 9.81 # m/s^2
    
    # Geometry calculations
    A = W * D # Cross-sectional area
    P = 2 * (W + D) # Perimeter
    Dh = 4 * A / P # Hydraulic diameter
    
    # Air properties at T_inf (20C = 293.15K) approx
    # Let's refine properties using film/bulk temp or just standard air approx.
    # To be accurate, we iterate. Let's write an iteration loop.
    # We need to find V_avg. 
    # Guess V_avg, find T_out, find T_avg, re-evaluate properties.
    
    # Properties function for air as a function of T (in K)
    def air_properties(T):
        # T in Kelvin
        # Simple polynomial fits or standard linear approximations around 20-30C
        # T_ref = 293.15 (20C)
        # rho = p / (R*T), p = 101325, R = 287.05
        rho = 101325 / (287.05 * T)
        # cp approx 1005 J/kgK
        cp = 1005.0
        # Dynamic viscosity mu (Sutherland's law or linear)
        mu = 1.81e-5 * (T / 293.15)**0.7
        nu = mu / rho
        # Thermal conductivity k
        k = 0.0257 * (T / 293.15)**0.8
        alpha = k / (rho * cp)
        Pr = nu / alpha
        beta = 1.0 / T # Ideal gas
        return rho, cp, mu, nu, k, alpha, Pr, beta

    # Iteration loop
    V_avg = 0.5 # initial guess m/s
    rho_inf, cp_inf, _, _, _, _, _, _ = air_properties(T_inf)
    
    for iteration in range(100):
        # 1. Estimate fluid properties at mean bulk/film temperature or iterate T_out
        # Let's assume an initial bulk mean temperature to compute heat transfer
        # We need the Nusselt number. 
        # For a vertical plate or channel, Ra_H is based on H.
        # Ra_H = g * beta * (Tw - T_inf) * H^3 / (nu * alpha)
        # Let's evaluate properties at film temperature T_film = (Tw + T_inf) / 2
        T_film = (Tw + T_inf) / 2.0
        rho_f, cp_f, mu_f, nu_f, k_f, alpha_f, Pr_f, beta_f = air_properties(T_film)
        
        Ra_H = (g * beta_f * (Tw - T_inf) * H**3) / (nu_f * alpha_f)
        
        # Churchill and Chu correlation for Nusselt number on a vertical plate
        Nu_H = (0.825 + (0.387 * Ra_H**(1/6)) / (1 + (0.492 / Pr_f)**(9/16))**(8/27))**2
        h_bar = Nu_H * k_f / H
        
        # Total heat transfer area = P * H
        As = P * H
        
        # Now, with V_avg, mass flow rate is m_dot = rho * A * V_avg
        # But rho varies. Let's use rho_inf for the inlet mass flow rate balance.
        m_dot = rho_inf * A * V_avg
        
        # Heat transferred to the fluid: Q_dot = h_bar * As * Delta_T_lm or simple energy balance
        # For a heated channel with constant wall temperature:
        # T_out = Tw - (Tw - T_inf) * exp(- (h_bar * As) / (m_dot * cp_f))
        # Let's check if the exponent is valid
        if m_dot > 0:
            exponent = (h_bar * As) / (m_dot * cp_f)
            T_out = Tw - (Tw - T_inf) * np.exp(-exponent)
        else:
            T_out = Tw
            
        T_mean = (T_inf + T_out) / 2.0
        rho_mean, _, _, nu_mean, _, _, _, _ = air_properties(T_mean)
        
        # Buoyancy driving pressure: Delta_p = (rho_inf - rho_mean) * g * H
        Delta_p_buoyancy = (rho_inf - rho_mean) * g * H
        
        # Friction factor (turbulent or laminar)
        Re_Dh = (V_avg * Dh) / nu_mean
        
        if Re_Dh < 2300:
            f = 64 / Re_Dh # basic laminar friction factor (can be adjusted for aspect ratio, but let's see if it's turbulent)
        else:
            # Swamee-Jain formula for smooth pipes (epsilon = 0)
            f = 0.25 / (np.log10(5.74 / Re_Dh**0.9))**2
            
        K_in = 0.5
        K_out = 1.0
        
        # Pressure drop: Delta_p_drop = [K_in + K_out + f * (H / Dh)] * (rho_mean * V_avg^2 / 2)
        # Balance: Delta_p_buoyancy = Delta_p_drop
        # Solve for new V_avg
        denom = (K_in + K_out + f * (H / Dh)) * (rho_mean / 2.0)
        if denom > 0 and Delta_p_buoyancy > 0:
            V_avg_new = np.sqrt(Delta_p_buoyancy / denom)
        else:
            V_avg_new = 0.001
            
        if abs(V_avg_new - V_avg) < 1e-5:
            V_avg = V_avg_new
            break
        V_avg = 0.5 * V_avg + 0.5 * V_avg_new

    # Final metrics
    Q_heat = m_dot * cp_f * (T_out - T_inf) # Watts
    # Enthalpy flow relative to ambient temperature reference (h = cp * (T - T_inf))
    enthalpy_flow = m_dot * cp_f * (T_out - T_inf) 
    
    fds_dir = Path(__file__).resolve().parents[3]
    results_dir = fds_dir.parent / 'out' / 'Convection'
    fds_results = []
    for resolution in ('5mm', '1cm', '2cm'):
        infile = results_dir / f'hot_chimney_{resolution}_devc.csv'
        with infile.open(newline='') as fid:
            next(fid)  # Units; device names are on the second row.
            reader = csv.DictReader(fid)
            last_row = None
            for last_row in reader:
                pass
        if last_row is None:
            raise ValueError(f'No device data in {infile}')
        fds_results.append(last_row)

    # Use the original printed precision; convert correlation enthalpy to kW.
    rows = [
        (r'Heat Transfer Coefficient ($\overline{h}$, \unit{W/(m^2 K)})', 'h_avg', h_bar, 2),
        (r'Average Velocity ($V_{\rm avg}$, m/s)', 'w_avg', V_avg, 4),
        (r'Air Exit Temperature ($T_{\rm out}$, \unit{\degreeCelsius})', 'T_avg', T_out - 273.15, 2),
        (r'Mass Flow Rate ($\dot{m}$, kg/s)', 'm_dot', m_dot, 4),
        (r'Total Enthalpy Flow Rate ($\dot{Q}$, kW)', 'q_dot', Q_heat / 1000.0, 2),
    ]
    table = [
        r'\begin{table}[!h]',
        r'\centering',
        r'\caption[Results of the \ct{hot_chimney} test cases]{Results of the \ct{hot_chimney} test cases.}',
        r'\label{hot_chimney_table}',
        r'\begin{tabular}{l|c||ccc}',
        r'\hline',
        r'Quantity & Correlation & FDS (5 mm) & FDS (1 cm) & FDS (2 cm) \\ \hline',
    ]
    for label, device, correlation, precision in rows:
        values = [correlation] + [float(result[device]) for result in fds_results]
        table.append(label + ' & ' + ' & '.join(f'{value:.{precision}f}' for value in values) + r' \\')
    table.extend([r'\hline', r'\end{tabular}', r'\end{table}'])

    outfile = fds_dir / 'Manuals' / 'FDS_Verification_Guide' / 'hot_chimney.tex'
    outfile.write_text('\n'.join(table) + '\n', encoding='utf-8')

hot_chimney()

