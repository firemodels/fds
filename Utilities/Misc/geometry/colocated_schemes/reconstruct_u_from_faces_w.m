function u_cell = reconstruct_u_from_faces_w(Uface, cells, faces, wtype, beta)
if nargin<4, wtype='trans'; end
if nargin<5, beta=0.3; end

Nc = numel(cells); dim = numel(cells(1).xc);
u_cell = zeros(Nc, dim);

for p = 1:Nc
    M = zeros(dim,dim); b = zeros(dim,1);
    flist = cells(p).faces(:);
    for k = 1:numel(flist)
        f  = flist(k);
        nf = faces(f).nf(:);      % owner→neighbor
        Af = faces(f).Af;

        % outward normal for this cell:
        s_pf = +1; if faces(f).neigh == p, s_pf = -1; end
        nout = s_pf * nf;

        % choose weight
        if isfield(faces,'dPN') && ~isempty(faces(f).dPN)
            dPN = faces(f).dPN;
        else
            % fallback estimate (only used for weighting)
            oth = faces(f).owner + faces(f).neigh - p;
            dvec = cells(oth).xc(:) - cells(p).xc(:);
            dPN  = abs(dvec.'*nout);
        end
        switch lower(wtype)
            case 'area'    , wf = Af;
            case 'trans'   , wf = Af/max(dPN,eps);
            case 'diamond' , wf = Af*dPN;
            case 'blend'   , wf = (1-beta)*Af + beta*Af/max(dPN,eps);
            otherwise      , wf = Af;
        end

        % --- Correct RHS build (key fix) ---
        % EITHER:   b = b + wf * (s_pf*Uface(f)) * nout;
        % OR (equivalent, safer): 
        b = b + wf * Uface(f) * nf;     % uses nf, not nout

        M = M + wf * (nout * nout.');
    end
    trM = trace(M); if trM<=0, trM=1; end
    M = M + (1e-12*trM)*eye(dim);
    u_cell(p,:) = (M\b).';
end
end
