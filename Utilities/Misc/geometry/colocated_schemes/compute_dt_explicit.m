function out = compute_dt_explicit(cells, faces, Uface, mu, rho, CFL, VN, Lx, Ly)
%COMPUTE_DT_EXPLICIT  Stable dt for Explicit Euler on unstructured FV mesh.
%
% dt_adv_p = CFL * Vp / sum_f |U_f| Af
% dt_dif_p = VN  * Vp / sum_f ( 2 * nu_f * Af / dn )
%
% Inputs:
%   cells : struct array with fields .V, .xc(1x2), .faces(1:nfp)
%   faces : struct array with fields .owner, .neigh, .nf(1x2 unit), .Af, .cf(1x2)
%   Uface : Nf x 1 face-normal velocities (along faces(f).nf, owner->neigh)
%   mu    : Nc x 1 dynamic viscosity at cells
%   rho   : Nc x 1 density at cells
%   CFL   : scalar convective Courant number (e.g. 0.5 … 0.9)
%   VN    : scalar Von Neumann number for diffusion (e.g. 0.5 is safe)
%   Lx,Ly : domain lengths (for periodic minimum-image distances)
%
% Output 'out' struct:
%   out.dt_adv_min   : min convective dt over cells
%   out.dt_dif_min   : min diffusive  dt over cells
%   out.dt           : min(out.dt_adv_min, out.dt_dif_min)
%   out.dt_adv_cell  : Nc x 1 per-cell convective dt
%   out.dt_dif_cell  : Nc x 1 per-cell diffusive  dt
%   out.limiter_cell : Nc x 1: -1 advective-limited, +1 diffusive-limited, 0 neither
%
% Notes:
% - Assumes every face has owner & neigh (periodic). For boundary faces w/ neigh=0,
%   treat dn and face-averages accordingly (not needed for fully periodic meshes).
% - Make sure Uface uses the same face-normal orientation as faces(f).nf.

Nc = numel(cells);
Nf = numel(faces);

% Precompute per-cell accumulators
sum_absU_A = zeros(Nc,1);      % Σ |U_f| Af
sum_dif    = zeros(Nc,1);      % Σ 2 * nu_f * Af / dn

% Handy inline for periodic minimum-image delta
wrap = @(dr,L)(dr - L.*(dr > 0.5*L) + L.*(dr <= -0.5*L));

for f = 1:Nf
    P  = faces(f).owner;
    N  = faces(f).neigh;                  % periodic mesh => N>0
    nf = faces(f).nf(:);
    Af = faces(f).Af;

    % Convective accumulator (same for both adjacent cells, with opposite sign handled later):
    absU_A = abs(Uface(f)) * Af;

    % Diffusion: face nu and compact dn
    if N > 0
        % arithmetic face average (ok for ν)
        nu_f = 0.5*(mu(P)/rho(P) + mu(N)/rho(N));
        dvec = [cells(N).xc(1) - cells(P).xc(1),  cells(N).xc(2) - cells(P).xc(2)];
        dvec = [wrap(dvec(1),Lx), wrap(dvec(2),Ly)];
        dn   = max( eps, dot(dvec, nf) );         % ensure positive
    else
        % boundary: one-sided estimate (not used for periodic)
        nu_f =  (mu(P)/rho(P));
        % estimate dn from owner-centroid to face-centroid projected on nf
        dvec = [faces(f).cf(1) - cells(P).xc(1), faces(f).cf(2) - cells(P).xc(2)];
        dn   = max( eps, dot(dvec, nf) );
    end

    dif_f = 2.0 * nu_f * (Af / dn);

    % Add to owner and neighbor cells
    sum_absU_A(P) = sum_absU_A(P) + absU_A;
    sum_dif(P)    = sum_dif(P)    + dif_f;

    if N > 0
        sum_absU_A(N) = sum_absU_A(N) + absU_A;
        sum_dif(N)    = sum_dif(N)    + dif_f;
    end
end

Vp = reshape([cells.V],[],1);

% Per-cell dt’s (protect against division by zero)
dt_adv_cell = inf(Nc,1);
dt_dif_cell = inf(Nc,1);

mask_adv = sum_absU_A > 0;
dt_adv_cell(mask_adv) = CFL * Vp(mask_adv) ./ sum_absU_A(mask_adv);

mask_dif = sum_dif > 0;
dt_dif_cell(mask_dif) = VN  * Vp(mask_dif) ./ sum_dif(mask_dif);

% Outputs
out.dt_adv_cell  = dt_adv_cell;
out.dt_dif_cell  = dt_dif_cell;
out.dt_adv_min   = min(dt_adv_cell);
out.dt_dif_min   = min(dt_dif_cell);
out.dt           = min(out.dt_adv_min, out.dt_dif_min);

% Limiter label per cell
out.limiter_cell = zeros(Nc,1);
adv_better = dt_adv_cell < dt_dif_cell;
dif_better = dt_dif_cell < dt_adv_cell;
out.limiter_cell(adv_better) = -1;   % advective-limited
out.limiter_cell(dif_better) = +1;   % diffusive-limited
end
