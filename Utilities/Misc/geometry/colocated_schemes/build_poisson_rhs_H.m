function rhs = build_poisson_rhs_H(U_f, divU_th, dt, cells, faces, FA_face_n, FB_face_n, Fextra_face_n)
% U_f         : Nf x 1  face-normal velocity (owner→neighbor)
% divU_th     : Nc x 1  target thermodynamic divergence per cell
% dt          : scalar
% FA_face_n   : Nf x 1  face-normal advective+body (scalar). [] → zeros
% FB_face_n   : Nf x 1  face-normal -p ∂n(1/ρ). [] → zeros
% Fextra_face_n: Nf x 1 optional extra face term (e.g., explicit viscous). [] → zeros
% Output:
% rhs(p) = - [ (div_th*V - Σ U_f A_f)/dt + Σ FA_face_n A_f + Σ FB_face_n A_f + Σ Fextra A_f ]

Nc = numel(cells);  Nf = numel(faces);
if nargin<6 || isempty(FA_face_n),    FA_face_n    = zeros(Nf,1); end
if nargin<7 || isempty(FB_face_n),    FB_face_n    = zeros(Nf,1); end
if nargin<8 || isempty(Fextra_face_n),Fextra_face_n= zeros(Nf,1); end

rhs = zeros(Nc,1);

for p = 1:Nc
    Vp = cells(p).V;
    S_U  = 0.0;    % Σ s_pf U_f A_f
    S_FA = 0.0;    % Σ s_pf F_A,f A_f
    S_FB = 0.0;    % Σ s_pf F_B,f A_f
    S_FX = 0.0;    % Σ s_pf Fextra,f A_f

    fids = cells(p).faces;
    for k = 1:numel(fids)
        f  = fids(k);
        Af = faces(f).Af;

        s = +1; if faces(f).neigh==p, s = -1; end  % outward sign for cell p

        S_U  = S_U  + s * U_f(f)        * Af;
        S_FA = S_FA + s * FA_face_n(f)  * Af;
        S_FB = S_FB + s * FB_face_n(f)  * Af;
        S_FX = S_FX + s * Fextra_face_n(f) * Af;
    end

    rhs(p) = - ( (divU_th(p)*Vp - S_U)/dt + S_FA + S_FB + S_FX );
end
end
