function [FA_cell, FA_face_n] = build_adv_body_forces_2d(u_in, rho, rho0, gvec, U_f, cells, faces, stage)
% F_A = div(u⊗u) - grad(0.5|u|^2) - (div u) u  - (1/rho)(rho-rho0) g
% u   : Nc x 2  (cell-centered)
% rho : Nc x 1
% U_f : Nf x 1  (face-normal velocity along faces(f).nf, owner->neigh)
% returns:
%   FA_cell   : Nc x 2  (cell-centered vector)
%   FA_face_n : Nf x 1  (face-normal avg of F_A for Poisson RHS)

global H_projection

%if stage==2
%   u = reconstruct_u_from_faces_w(U_f, cells, faces);
%else
   u = u_in;
%end
Nc = size(u,1); Nf = numel(faces);
Kc = 0.5*sum(u.^2,2);           % 0.5|u|^2 at cells



sum_divuV = zeros(Nc,1);        % Σ s_pf U_f A_f
sum_divUU = zeros(Nc,2);        % Σ s_pf (u_f^c U_f) A_f
sum_gradKV= zeros(Nc,2);        % Σ s_pf (K_f^c n_f) A_f

for f = 1:Nf
    P = faces(f).owner; N = faces(f).neigh;
    nf = faces(f).nf(:); Af = faces(f).Af;
    Uf = U_f(f);

    % centered face values (periodic mesh ⇒ N>0)
    if N>0
        uf = 0.5*(u(P,:)+u(N,:));
        Kf = 0.5*(Kc(P)+Kc(N));
    else
        uf = u(P,:);  Kf = Kc(P);
    end

    adv_flux   = (uf .* Uf) * Af;            % vector
    gradK_flux = (Kf * nf.') * Af;           % vector

    % owner (+), neighbor (-)
    sum_divUU(P,:) = sum_divUU(P,:) + adv_flux;
    sum_gradKV(P,:)= sum_gradKV(P,:) + gradK_flux;
    sum_divuV(P)   = sum_divuV(P)   + Uf*Af;

    if N>0
        sum_divUU(N,:) = sum_divUU(N,:) - adv_flux;
        sum_gradKV(N,:)= sum_gradKV(N,:) - gradK_flux;
        sum_divuV(N)   = sum_divuV(N)   - Uf*Af;
    end
end

FA_cell = zeros(Nc,2);
g = gvec(:).';
for p = 1:Nc
    Vp     = cells(p).V;
    divUU  = sum_divUU(p,:)/Vp;              % ∇·(u⊗u)
    if(H_projection)
       gradK  = sum_gradKV(p,:)/Vp;          % ∇(0.5|u|^2)
       divu_u = (sum_divuV(p)/Vp)*u(p,:);    % (∇·u) u
    else
       gradK  = zeros(1,2);
       divu_u = zeros(1,2);
    end
    body   = - ((rho(p)-rho0)/rho(p)) * g;
    FA_cell(p,:) = divUU - gradK - divu_u + body;
end

% face-normal average of F_A (for Poisson RHS)
FA_face_n = zeros(Nf,1);
for f = 1:Nf
    P = faces(f).owner; N = faces(f).neigh; nf = faces(f).nf(:);
    Fav = (N>0) * 0.5*(FA_cell(P,:)+FA_cell(N,:)) + (N==0)*FA_cell(P,:);
    FA_face_n(f) = Fav * nf;
end
end
