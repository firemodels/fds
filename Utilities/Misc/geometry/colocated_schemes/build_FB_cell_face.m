function [FB_cell, FB_face_n] = build_FB_cell_face(p_cell, rho_cell, cells, faces, Lx, Ly)
% p_cell   : Nc x 1   pressure (or Bernoulli-pressure part) at cells
% rho_cell : Nc x 1   density at cells
% Outputs:
%   FB_cell   : Nc x 2   vector  -p * ∇(1/ρ)  at cells
%   FB_face_n : Nf x 1   scalar  face-normal avg: 0.5(FB_P+FB_N)·n_f

Nc = numel(cells); Nf = numel(faces);

invR = 1 ./ rho_cell;                     % 1/ρ at cells
gradInvR = ls_grad_scalar_2d(invR, cells, faces, Lx, Ly);   % ∇(1/ρ) at cells

% Cell-centered F_B = - p * ∇(1/ρ)
FB_cell = - p_cell(:) .* gradInvR;        % implicit expansion: (Nc×1) .* (Nc×2) → (Nc×2)

% Face-normal average like F_A
FB_face_n = zeros(Nf,1);
for f = 1:Nf
    P  = faces(f).owner; 
    N  = faces(f).neigh;                  % periodic mesh ⇒ N>0
    nf = faces(f).nf(:);
    if N>0
        FBav = 0.5*(FB_cell(P,:) + FB_cell(N,:));
    else
        FBav = FB_cell(P,:);              % (kept for completeness)
    end
    FB_face_n(f) = FBav * nf;
end
end

% function FB_face_n = build_FB_face_n_compact(p, rho, cells, faces, Lx, Ly)
% Nf = numel(faces);
% FB_face_n = zeros(Nf,1);
% for f = 1:Nf
%     P  = faces(f).owner; 
%     N  = faces(f).neigh;              % periodic ⇒ N>0
%     nf = faces(f).nf(:);
% 
%     dvec = periodic_delta_2d(cells(N).xc - cells(P).xc, Lx, Ly);
%     dn   = dot(dvec, nf);             % compact face distance along nf
%     invrP = 1.0 / rho(P);
%     invrN = 1.0 / rho(N);
%     d_invrho_dn = (invrN - invrP) / dn;
% 
%     pf = 0.5*(p(P) + p(N));           % centered face pressure
%     FB_face_n(f) = - pf * d_invrho_dn;
% end
% end
