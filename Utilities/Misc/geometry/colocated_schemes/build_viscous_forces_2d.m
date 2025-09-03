function [divTau_cell, FA_visc_cell, FA_visc_face_n,traction_f] = build_viscous_forces_2d(u, rho, mu, lambda, cells, faces, Lx, Ly)
% Compute (∇·tau) at cells using face tractions of the Newtonian viscous stress.
% u:      Nc x 2  velocities (u,v) at cells
% rho:    Nc x 1   density
% mu:     Nc x 1  dynamic viscosity at cells (can vary)
% lambda: Nc x 1  bulk viscosity at cells (if [], uses Stokes: lambda = -2/3 mu)
% cells:  struct array with .xc(1x2), .V, .faces(1xnf)
% faces:  struct array with .owner, .neigh (0 if boundary), .nf(1x2 unit), .Af, .cf(1x2)
%
% returns:
% divTau_cell: Nc x 2  (∇·tau) at cells
% FA_visc_face_n : Nf x 1  face-normal scalar to add into FA_face_n
% traction_f:  Nf x 2  tau_f * n_f (vector without area)

Nc = size(u,1);
Nf = numel(faces);

if isempty(lambda)
    lambda = -2/3 * mu;   % Stokes' hypothesis
end

% 1) cell-centered velocity gradients by LS
[gradU, gradV] = ls_grad_velocity_2d(u, cells, faces, Lx, Ly);
% gradU(p,:) = [du/dx, du/dy]; gradV(p,:) = [dv/dx, dv/dy]

% 2) face traction tau_f * n_f
traction_f = zeros(Nf,2);

for f = 1:Nf
    P  = faces(f).owner;
    N  = faces(f).neigh;
    nf = faces(f).nf(:);      % unit normal (P -> N)
    % build face gradient by centered average (robust, symmetric)
    if N>0
        Gp = [gradU(P,1), gradU(P,2);   % ∇u at P (2x2)
              gradV(P,1), gradV(P,2)];  % rows: du/dx du/dy ; dv/dx dv/dy
        Gn = [gradU(N,1), gradU(N,2);
              gradV(N,1), gradV(N,2)];
        Gf = 0.5*(Gp + Gn);

        muf = 0.5*(mu(P) + mu(N));
        lamf= 0.5*(lambda(P) + lambda(N));
    else
        % boundary face: simple one-sided fallback using cell gradient
        Gf  = [gradU(P,1), gradU(P,2);
               gradV(P,1), gradV(P,2)];
        muf = mu(P);
        lamf= lambda(P);
        % NOTE: For wall/inlet/outlet you can/should replace this
        % with BC-consistent gradient (ghost-cell or analytic).
    end

    % divergence at face (trace of gradient)
    divUf = Gf(1,1) + Gf(2,2);

    % stress tensor at face: tau = mu*(Grad u + Grad u^T) + lambda*(div u) I
    tau = muf*(Gf + Gf.') + lamf*divUf*eye(2);

    % traction: t_f = tau * n_f   (vector 2x1)
    traction_f(f,:) = (tau * nf).';
end

% 3) assemble (∇·tau)_cell by divergence theorem
divTau_cell = zeros(Nc,2);
for p = 1:Nc
    sumF = [0,0];
    fids = cells(p).faces;
    for k = 1:numel(fids)
        f  = fids(k);
        Af = faces(f).Af;
        nf = faces(f).nf(:);    % orientation P->N

        if faces(f).owner == p
            % outward normal of cell p is +nf
            sumF = sumF + (traction_f(f,:)*Af);
        elseif faces(f).neigh == p
            % outward normal of cell p is -nf
            sumF = sumF - (traction_f(f,:)*Af);
        else
            % should not happen
        end
    end
    divTau_cell(p,:) = sumF / cells(p).V;
end

% Cell viscous part of F_A (vector)
FA_visc_cell = - divTau_cell ./ rho;   % Nc x 2, implicit expansion

% Face-normal average
FA_visc_face_n = zeros(Nf,1);
for f = 1:Nf
    P  = faces(f).owner;
    N  = faces(f).neigh;
    nf = faces(f).nf(:);
    if N>0
        Fav = 0.5 * (FA_visc_cell(P,:) + FA_visc_cell(N,:));
    else
        Fav = FA_visc_cell(P,:);
    end
    FA_visc_face_n(f) = Fav * nf;   % dot product
end
end