function [u_out, Uface_out, H_cell, p_out, dHdn_f, gradH] = ...
    do_projection_stage(u_in, u_in1, Uface_in, Uface_in1, p_in, rho, mu, lambda, ...
                        cells, faces, Lx, Ly, dt, rho0, gvec, divU_th, ...
                        A, poisson_opts, add_meanK, stage)

global H_projection

% ---------- build viscous ----------
[~, FA_visc_cell, FA_visc_face_n,~] = ...
    build_viscous_forces_2d(u_in1, rho, mu, lambda, cells, faces, Lx, Ly);

% ---------- advective + body (uses GIVEN face U at this stage) ----------
[FAdv_cell, FAdv_face_n] = ...
    build_adv_body_forces_2d(u_in1, rho, rho0, gvec, Uface_in, cells, faces, stage);

% ---------- baroclinic/B-term in cells then faces ----------
[FB_cell, FB_face_n] = build_FB_cell_face(p_in, rho, cells, faces, Lx, Ly);

% Combine face & cell contributions for RHS and cell update
FA_face_n = FAdv_face_n + FA_visc_face_n;  % faces
FA_cell   = FAdv_cell   + FA_visc_cell;    % cells

% ---------- Poisson RHS (flux form) ----------
rhs = build_poisson_rhs_H(Uface_in, divU_th, dt, cells, faces, FA_face_n, FB_face_n);

% ---------- Solve Poisson for H (your sign convention) ----------
[H_cell, dHdn_f, gradH] = poisson_solve(A, rhs, cells, faces, Lx, Ly, poisson_opts);

% ---------- (optional) gauge tweak to help p match analytic TG) ----------
if(H_projection)
    if add_meanK
        K_cell = 0.5*sum(u_in.^2,2);
        H_gauge = H_cell + mean(K_cell);
    else
        H_gauge = H_cell;
    end
end

% ---------- Project faces (mass consistency) ----------
Uface_out = project_faces(Uface_in, Uface_in1, dt, FA_face_n, FB_face_n, dHdn_f, stage);

% ---------- Project cells (momentum update) ----------
if(stage==1)
   u_out = reconstruct_u_from_faces_w(Uface_out, cells, faces);
%   u_out2 = project_cells(u_in, dt, FA_cell, divTau_cell, rho, FB_cell, gradH);
%   u_out  = 0.5*(u_out1 + u_out2);
else
    u_out = project_cells(u_in, u_in1, dt, FA_cell, FB_cell, gradH, stage);
end

if(H_projection)
    % ---------- Pressure from Bernoulli for NEXT stage inputs ----------
    K_out = 0.5*sum(u_out.^2,2);
    p_out = rho .* (H_gauge - K_out);
else
    p_out = H_cell;
end

end
