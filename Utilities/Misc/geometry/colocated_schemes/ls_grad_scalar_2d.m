function gradPhi = ls_grad_scalar_2d(phi, cells, faces, Lx, Ly)
% phi: Nc x 1 cell scalar
Nc = numel(cells);
gradPhi = zeros(Nc,2);
for p = 1:Nc
    xp  = cells(p).xc;
    fids = cells(p).faces;
    M = zeros(2,2); b = [0;0];
    for k = 1:numel(fids)
        f = fids(k);
        P = faces(f).owner; N = faces(f).neigh;
        if     P==p && N>0, q = N;
        elseif N==p && P>0, q = P;
        else, continue; end   % (no boundary faces in periodic mesh)
        rq  = periodic_delta_2d(cells(q).xc - xp, Lx, Ly);
        Af  = faces(f).Af;
        w   = Af / (dot(rq,rq) + 1e-14);     % robust weight
        M = M + w * (rq(:)*rq(:).');
        b = b + w * (phi(q) - phi(p)) * rq(:);
    end
    trM = trace(M); if trM<=0, trM = 1; end
    M = M + (1e-12*trM)*eye(2);             % tiny regularization
    gradPhi(p,:) = (M\b).';
end
end
