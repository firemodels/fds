function [u, Uface, p, H_cell, t, hist] = run_ssprk2( ...
    u, Uface, p, rho, mu, lambda, ...
    cells, faces, Lx, Ly, params)
global compact_laplacian Nx Ny Xc Yc dstep_plot do_ssprk2 H_projection
% ----------------- defaults -----------------
Nc = numel(cells);
defaults = struct( ...
    'rho0',      1.0, ...
    'gvec',      [0,0], ...
    'CFL',       0.8, ...
    'VN',        0.5, ...
    'divU_th',   zeros(Nc,1), ...
    't_begin',   0.0, ...
    't_end',     1.0, ...
    'dt',        0.01, ...
    'variable_time_step', 1, ...
    'max_steps', 10^6, ...
    'add_meanK', true, ...
    'poisson_opts', struct('fix','pin','pin_idx',1,'use_ilu',true), ...
    'tg_params', [] ...    % struct with fields U0,kx,ky,nu (optional)
    );
params = set_defaults(params, defaults);

% Prebuild compact Poisson matrix
A = build_poisson_matrix_compact(cells, faces, Lx, Ly);

% History
hist = struct('step', [], 't',[], 'dt',[], 'L2_mass',[], ...
              'L2_m_abs',[], 'L2_m_rel',[], 'M_num',[], 'M_ex',[], ...
              'L2_k_abs',[], 'L2_k_rel',[], 'K_num_tot',[], ...
              'K_ex_tot',[], 'K_ex_formula',[], 'u_num_mid', [], ...
              'u_ex_mid', []);

% Initial time parameters:
t = params.t_begin;
dt= params.dt;
step = 0;

% Initial Mass, Momentum and TKE:
divU_after = cell_divergence_from_faces(Uface, cells, faces);
res_mass   = divU_after - params.divU_th;
L2_mass    = norm(res_mass,2);
tg = params.tg_params;  % struct U0,kx,ky,nu
% Momentum
[L2m_abs, L2m_rel, M_num, M_ex] = l2_error_momentum(u, rho, cells, t, tg);
% KE
[L2k_abs, L2k_rel, K_num_tot, K_ex_tot, K_ex_formula] = ...
            l2_error_ke(u, rho, cells, t, tg);

% ---- Save history ----
hist.step(end+1,1)    = step;
hist.t(end+1,1)       = t;
hist.dt(end+1,1)      = dt;
hist.L2_mass(end+1,1) = L2_mass;

hist.L2_m_abs(end+1,1)= L2m_abs;
hist.L2_m_rel(end+1,1)= L2m_rel;
hist.M_num(end+1,1:2) = M_num;
hist.M_ex(end+1,1:2)  = M_ex;

hist.L2_k_abs(end+1,1)= L2k_abs;
hist.L2_k_rel(end+1,1)= L2k_rel;
hist.K_num_tot(end+1,1)  = K_num_tot;
hist.K_ex_tot(end+1,1)   = K_ex_tot;
hist.K_ex_formula(end+1,1)= K_ex_formula;

hist.u_num_mid(end+1,1:2)=u(params.p_mid,1:2);

% Optional short log
fprintf(['SSPRK2 %5d  t=%.6e  dt=%.3e  ||mass||_2=%.3e  ', ...
         'L2m=%.2e  L2k=%.2e  Knum=%.3e  Kex=%.3e\n'], ...
          step, t, dt, L2_mass, L2m_abs, L2k_abs, K_num_tot/(Lx*Ly), K_ex_tot/(Lx*Ly));

[fid]=plot_TG_field(t,Nx,Ny,Xc,Yc,u,p);
figure(gcf)
pause(0.2)

% Time loop
while (t < params.t_end) && (step < params.max_steps)
    step = step + 1;

    % ---- dt from CFL/VN ----
    if (params.variable_time_step==1) 
    dts = compute_dt_explicit(cells, faces, Uface, mu, rho, params.CFL, params.VN, Lx, Ly);
    else 
    dts.dt = dt;
    end
    dt  = min(dts.dt, params.t_end - t);

    % ---- Stage 1 ----
    [u1, Uface1, H1, p1] = do_projection_stage(u,u,Uface,Uface,p, rho, mu, lambda, ...
                         cells, faces, Lx, Ly, dt, ...
                         params.rho0, params.gvec, params.divU_th, ...
                         A, params.poisson_opts, params.add_meanK,1);

    if(do_ssprk2==1) % SSPRK2
    % ---- Stage 2 ----
    [u2, Uface2, H2, p2] = do_projection_stage(u,u1, ...
                         Uface,Uface1,p1, rho, mu, lambda, ...
                         cells, faces, Lx, Ly, dt, ...
                         params.rho0, params.gvec, params.divU_th, ...
                         A, params.poisson_opts, params.add_meanK,2);
    else % Forward Euler.
        u2=u1; H2=H1; Uface2=Uface1;
    end

    % ---- SSPRK2 combo ----
    u_next     = u2; 
    Uface_next = Uface2; 
    H_cell     = H2;  

    % Pressure from Bernoulli (with optional gauge tweak)
    if(H_projection)
        K_next = 0.5*sum(u_next.^2,2);
        if params.add_meanK
            H_eff = H_cell + mean(K_next);
        else
            H_eff = H_cell;
        end
        p_next = rho .* (H_eff - K_next);
    else
        H_eff  = H_cell;
        p_next = H_cell;
    end

    % ---- Mass residual ----
    divU_after = cell_divergence_from_faces(Uface_next, cells, faces);
    res_mass   = divU_after - params.divU_th;
    L2_mass    = norm(res_mass,2);

    % ---- Momentum / KE checks (if TG params available) ----
    if ~isempty(params.tg_params)
        tg = params.tg_params;  % struct U0,kx,ky,nu

        % Momentum
        [L2m_abs, L2m_rel, M_num, M_ex] = l2_error_momentum(u_next, rho, cells, t+dt, tg);

        % KE
        [L2k_abs, L2k_rel, K_num_tot, K_ex_tot, K_ex_formula] = ...
            l2_error_ke(u_next, rho, cells, t+dt, tg);
    else
        L2m_abs=NaN; L2m_rel=NaN; M_num=[NaN NaN]; M_ex=[NaN NaN];
        L2k_abs=NaN; L2k_rel=NaN; K_num_tot=NaN; K_ex_tot=NaN; K_ex_formula=NaN;
    end

    % ---- Commit & time advance ----
    u     = u_next;
    Uface = Uface_next;
    p     = p_next;
    t     = t + dt;

    % ---- check pressure work ------
    % check_product_rule_cellwise(u, H_eff, cells, faces)
    % check_H_divU(H_cell, Uface, cells, faces)

    % ---- Save history ----
    hist.step(end+1,1)    = step;
    hist.t(end+1,1)       = t;
    hist.dt(end+1,1)      = dt;
    hist.L2_mass(end+1,1) = L2_mass;

    hist.L2_m_abs(end+1,1)= L2m_abs;
    hist.L2_m_rel(end+1,1)= L2m_rel;
    hist.M_num(end+1,1:2) = M_num;
    hist.M_ex(end+1,1:2)  = M_ex;

    hist.L2_k_abs(end+1,1)= L2k_abs;
    hist.L2_k_rel(end+1,1)= L2k_rel;
    hist.K_num_tot(end+1,1)  = K_num_tot;
    hist.K_ex_tot(end+1,1)   = K_ex_tot;
    hist.K_ex_formula(end+1,1)= K_ex_formula;

    hist.u_num_mid(end+1,1:2)=u(params.p_mid,1:2);

    % Optional short log
    fprintf(['SSPRK2 %5d  t=%.6e  dt=%.3e  ||mass||_2=%.3e  ', ...
             'L2m=%.2e  L2k=%.2e  Knum=%.3e  Kex=%.3e\n'], ...
             step, t, dt, L2_mass, L2m_abs, L2k_abs, K_num_tot/(Lx*Ly), K_ex_tot/(Lx*Ly));

    % Plot tg field:
    conditional = mod(step,dstep_plot)<1e-3;
    if mod(step,dstep_plot)<1e-3
        [fid]=plot_TG_field(t,Nx,Ny,Xc,Yc,u,p);
        figure(gcf)
        pause(0.2)
    end
end
end

function s = set_defaults(s, d)
fn = fieldnames(d);
for k=1:numel(fn)
    if ~isfield(s,fn{k}) || isempty(s.(fn{k}))
        s.(fn{k}) = d.(fn{k});
    end
end
end

