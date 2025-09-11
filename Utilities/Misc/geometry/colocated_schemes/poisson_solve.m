function [H_cell, dHdn_f, gradH, info] = poisson_solve(A, rhs, cells, faces, Lx, Ly, opts)
%POISSON_SOLVE  Solve compact Poisson for H and return zero-mean H, face normals and LS gradients.
%
% Inputs
%   A      : Nc x Nc sparse Poisson matrix (compact Laplacian)
%   rhs    : Nc x 1 right-hand-side (flux form)
%   cells  : struct array with fields .V (scalar), .xc(1x2), .faces(1:nf)
%   faces  : struct array with fields .owner, .neigh, .nf(1x2 unit owner->neigh), .Af
%   Lx,Ly  : domain lengths (for periodic minimum-image deltas)
%   opts   : struct with optional fields:
%            .fix      = 'pin' (default) or 'mean'   % nullspace fix
%            .pin_idx  = 1 (default)                 % index to pin if fix='pin'
%            .tol      = 1e-12                       % PCG tolerance
%            .maxit    = 500                         % PCG max iterations
%            .use_ilu  = true                        % try ILU preconditioner
%            .ilu_droptol = 1e-6                     % ILU drop tolerance
%
% Outputs
%   H_cell : Nc x 1  solution with volume-weighted zero mean
%   dHdn_f : Nf x 1  compact face-normal derivative ( (H_N-H_P)/((xN-xP)·nf) )
%   gradH  : Nc x 2  LS reconstruction of ∇H from dHdn_f
%   info   : struct with .flag, .relres, .iter, .time_solve
%
% Notes
%   - Domain assumed fully periodic (pure Neumann); we remove the nullspace.
%   - For non-periodic cases, extend the matrix/RHS assembly accordingly.
global compact_laplacian

% Flag to set the gradient in cell centers either to divergence gradient (1) 
% or Least Squares(0):
do_divergence_gradient = 1;


Nc = numel(cells); Nf = numel(faces);

% ----- defaults
if nargin < 7, opts = struct(); end
opts = set_default(opts, 'fix',          'pin');
opts = set_default(opts, 'pin_idx',      1);
opts = set_default(opts, 'tol',          1e-12);
opts = set_default(opts, 'maxit',        5000);
opts = set_default(opts, 'use_ilu',      true);
opts = set_default(opts, 'ilu_droptol',  1e-6);

% ----- nullspace fix and change signs to make eigs of A > 0:
switch lower(opts.fix)
  case 'pin'
    [Afix, rhsfix] = pin_one_dof(-A, -rhs, opts.pin_idx);
    take = 1:Nc;
  case 'mean'
    [Afix, rhsfix] = mean_zero_constraint(-A, -rhs, cells);
    take = 1:Nc;   % first Nc entries are H; last one is Lagrange multiplier
  otherwise
    error('opts.fix must be ''pin'' or ''mean''.');
end

% ----- solve
t0 = tic;
L = []; U = [];
if opts.use_ilu
  try
    setup.type = 'ilutp';
    setup.droptol = opts.ilu_droptol;
    [L,U] = ilu(Afix, setup);
  catch
    L = []; U = [];
  end
end

try
  [Hsol, flag, relres, iter] = pcg(Afix, rhsfix, opts.tol, opts.maxit, L, U);
catch
  [Hsol, flag, relres, iter] = pcg(Afix, rhsfix, opts.tol, opts.maxit);
end
info.flag = flag; info.relres = relres; info.iter = iter; info.time_solve = toc(t0);
if flag ~= 0
  warning('poisson_solve: PCG did not fully converge (flag=%d, relres=%.2e, iter=%d).', flag, relres, iter);
end

% Extract H and enforce zero-mean
H_cell = Hsol(take);
V = reshape([cells.V], [], 1);
Hbar = sum(V .* H_cell) / sum(V);
H_cell = H_cell - Hbar;

if(compact_laplacian)

[dHdn_f, gradH] = grad_from_compact_normals_w(H_cell, cells, faces, Lx, Ly,'trans');
%[dHdn_f, gradH] = grad_from_compact_normals(H_cell, cells, faces, Lx, Ly);
% ----- face-normal compact derivative dH/dn
% dHdn_f = zeros(Nf,1);
% for f = 1:Nf
%   P  = faces(f).owner;
%   N  = faces(f).neigh;
%   nf = faces(f).nf(:);
%   dvec = periodic_delta_2d(cells(N).xc - cells(P).xc, Lx, Ly);
%   dn   = dot(dvec, nf);
%   dHdn_f(f) = (H_cell(N) - H_cell(P)) / dn;
% end
% 
% % ----- ∇H at cells -----
% gradH = zeros(Nc,2);
% 
% if do_divergence_gradient == 1
%     % === Green–Gauss (divergence theorem) using H_f = 0.5(H_P+H_N) ===
%     sum_flux = zeros(Nc,2);             % Σ H_f n_f A_f
%     for f = 1:Nf
%         P  = faces(f).owner;
%         N  = faces(f).neigh;
%         nf = faces(f).nf(:);
%         Af = faces(f).Af;
% 
%         % face-averaged H
%         if N > 0
%             Hf = 0.5*(H_cell(P) + H_cell(N));
%         else
%             % (boundary case; replace with BC value if needed)
%             Hf = H_cell(P);
%         end
% 
%         contrib = (Hf * nf.').*Af;      % 1x2 vector
% 
%         % owner (+), neighbor (−)
%         sum_flux(P,:) = sum_flux(P,:) + contrib;
%         if N > 0
%             sum_flux(N,:) = sum_flux(N,:) - contrib;
%         end
%     end
%     gradH = sum_flux ./ V;              % V is Nc×1 of cell volumes
% 
% else
%     % === Least–Squares from compact normals: M_p ∇H_p = Σ Af (∂H/∂n)_f n_f ===
%     for p = 1:Nc
%         M = zeros(2,2); b = [0;0];
%         fids = cells(p).faces;
%         for k = 1:numel(fids)
%             f   = fids(k);
%             nf  = faces(f).nf(:);
%             Af  = faces(f).Af;
%             qf  = dHdn_f(f);           % (H_N - H_P)/dn
%             M = M + Af * (nf*nf.');
%             b = b + Af * qf * nf;
%         end
%         trM = trace(M); if trM<=0, trM=1; end
%         M = M + (1e-12*trM)*eye(2);     % tiny Tikhonov for robustness
%         gradH(p,:) = (M\b).';
%     end
% end

else

    [gradH, ~, dHdn_f] = gradH_sparse_GG(H, cells, faces);

end

end

% ======= helpers =======

function [Afix, rhsfix] = pin_one_dof(A, rhs, pind)
Afix   = A;     rhsfix = rhs;
Afix(pind,:) = 0;  Afix(:,pind) = 0;
Afix(pind,pind) = 1;
rhsfix(pind) = 0;
end

function [Aaug, rhsaug] = mean_zero_constraint(A, rhs, cells)
Nc = numel(cells);
V  = reshape([cells.V], [], 1);
W  = V / sum(V);             % volume-weighted mean
Aaug  = [A,  W;  W.',  0];
rhsaug= [rhs; 0];
end

function s = set_default(s, name, val)
if ~isfield(s,name) || isempty(s.(name)), s.(name) = val; end
end
