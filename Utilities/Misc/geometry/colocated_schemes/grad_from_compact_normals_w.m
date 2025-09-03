function [dphidn_f, gradphi_cell] = grad_from_compact_normals_w(phi, cells, faces, Lx, Ly, wtype, beta)
if nargin<6 || isempty(wtype), wtype = 'trans'; end   % 'area' | 'trans' | 'diamond' | 'blend'
if nargin<7 || isempty(beta),  beta  = 0.3; end

Nf = numel(faces); Nc = numel(cells);
dphidn_f = zeros(Nf,1);


% face-normal derivatives from two-point differences
for f = 1:Nf
    P = faces(f).owner; N = faces(f).neigh; nf = faces(f).nf(:);
   
    if N>0
        dr_raw = cells(N).xc(:) - cells(P).xc(:);
        dr     = periodic_delta_2d(dr_raw, Lx, Ly);
        dphidn_f(f) = (phi(N) - phi(P)) / dot(dr,nf);
    else
        disp('N==0')
        dphidn_f(f) = 0.0; % set per BC if needed
    end
end

% LS gradient per cell with weights
gradphi_cell = zeros(Nc,2);
for p = 1:Nc
    M = zeros(2,2); b = zeros(2,1);
    flist = cells(p).faces(:);
    for k = 1:numel(flist)
        f  = flist(k);
        nf = faces(f).nf(:);
        Af = faces(f).Af;
        % outward normal for cell p
        s_pf = +1; if faces(f).neigh==p, s_pf = -1; end
        nout = s_pf*nf;

        % choose weight
        dr_raw = cells(N).xc(:) - cells(P).xc(:);
        dr     = periodic_delta_2d(dr_raw, Lx, Ly);
        dPN = abs(dot(dr,nf));
        switch lower(wtype)
            case 'area'    , wf = Af;
            case 'trans'   , wf = Af/max(dPN,eps);
            case 'diamond' , wf = Af*dPN;
            case 'blend'   , wf = (1-beta)*Af + beta*Af/max(dPN,eps);
            otherwise      , wf = Af;
        end

        qf = dphidn_f(f);      % owner->neigh normal derivative
        % Use owner->neigh normal in RHS (no extra sign):
        M = M + wf * (nout * nout.');
        b = b + wf * qf * nf;
    end
    %trM = trace(M); if trM<=0, trM=1; end
    %M = M + (1e-12*trM)*eye(2);
    gradphi_cell(p,:) = (M\b).';
end
end
