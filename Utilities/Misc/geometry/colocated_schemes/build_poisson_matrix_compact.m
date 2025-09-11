function A = build_poisson_matrix_compact(cells, faces, Lx, Ly)
% Returns sparse A (Nc x Nc) for: sum_f ((H_N - H_P)/dn)*Af = RHS

do_A_tests = 0; % If 1, do sanity tests on matrix.

Nc = numel(cells);
Nf = numel(faces);

% Compact (two-point) Poisson matrix in flux form
% A*phi produces per-cell SUM_f (Af/|dPN|) * (phi_P - phi_N), i.e., flux balance.
I = zeros(4*Nf,1); J = I; V = I; m = 0;

for f = 1:Nf
    P  = faces(f).owner;
    N  = faces(f).neigh;              % assume periodic/internal ⇒ N>0
    if N<=0, continue; end
    Af = faces(f).Af;
    nf = faces(f).nf(:);              % unit normal, owner→neighbor

    % owner→neighbor minimal-image displacement (periodic safe)
    dvec = periodic_delta_2d(cells(N).xc - cells(P).xc, Lx, Ly);
    dn   = dot(dvec, nf);             % can be ±
    %dPN  = max(eps, abs(dn));         % ALWAYS use absolute distance
    T    = Af / dn; %PN;                  % transmissibility (≥0)

    % symmetric 2x2 face stencil:
    % row P: +T at P, -T at N
    m=m+1; I(m)=P; J(m)=P; V(m)= T;
    m=m+1; I(m)=P; J(m)=N; V(m)=-T;
    % row N: -T at P, +T at N
    m=m+1; I(m)=N; J(m)=N; V(m)= T;
    m=m+1; I(m)=N; J(m)=P; V(m)=-T;
end

A = -sparse(I(1:m), J(1:m), V(1:m), Nc, Nc);


% Tests:
if(do_A_tests)
% symmetry
fprintf('||A - A''|| = %.3e\n', norm(A - A','fro'));

% nullspace (periodic): one zero eigenvalue (the constant)
lam = eigs(A,4,'smallestreal');
disp(sort(real(lam)).');   % expect ~[0, +, +, +]

z = A*ones(numel(cells),1);
fprintf('||A*1||_2 = %.3e\n', norm(z,2));


% divergence-of-gradient pairing (with compact normals)
H = rand(Nc,1); H = H - mean(H);
[dHdn_f, ~] = grad_from_compact_normals(H, cells, faces, Lx, Ly);
LHS = A*H;                      % face flux balance per cell
RHS = zeros(Nc,1);
for p=1:Nc
  acc=0; for f=cells(p).faces
    s = +1; if faces(f).neigh==p, s=-1; end
    acc = acc + s*dHdn_f(f) * faces(f).Af;
  end
  RHS(p) = acc;
end
fprintf('pairing ||A*H - sum_f (dHdn_f Af s_pf)||_2 = %.3e\n', norm(LHS-RHS,2));

bad = 0;
for f=1:Nf
  P=faces(f).owner; N=faces(f).neigh;
  dvec = periodic_delta_2d(cells(N).xc - cells(P).xc, Lx, Ly);
  if dot(dvec, faces(f).nf(:)) <= 0, bad=bad+1; end
end
fprintf('faces with nf·(xN-xP) <= 0: %d\n', bad);

% 1) normals point owner->neighbor on average (not required to be perfectly aligned)
bad = 0;
for f = 1:numel(faces)
    P=faces(f).owner; N=faces(f).neigh;
    if N<=0, continue; end
    dvec = periodic_delta_2d(cells(N).xc - cells(P).xc, Lx, Ly);
    if dot(dvec, faces(f).nf(:)) < 0, bad = bad + 1; end
end
fprintf('faces with (xN-xP)·nf < 0: %d\n', bad);
% It’s OK if some are negative (skew), but we rely on ABS(dPN) in both A and q_out.

% 2) row-sum (nullspace) check: should be 0 for unpinned periodic A
fprintf('||A*1||_2 = %.3e\n', norm(A*ones(numel(cells),1),2));
end
end
