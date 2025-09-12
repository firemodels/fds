function [gradU, gradV] = ls_grad_velocity_2d(u, cells, faces, Lx, Ly)
Nc = size(u,1);
gradU = zeros(Nc,2); gradV = zeros(Nc,2);

for p = 1:Nc
    xp = cells(p).xc;
    fids = cells(p).faces;
    M = zeros(2,2); bu = [0;0]; bv = [0;0];

    for k = 1:numel(fids)
        f = fids(k);
        P = faces(f).owner; N = faces(f).neigh;
        if     P==p && N>0, q = N;
        elseif N==p && P>0, q = P;
        else,               continue
        end
        rq_raw = cells(q).xc - xp;
        rq     = periodic_delta_2d(rq_raw, Lx, Ly);      % <-- wrap here
        Af = faces(f).Af;
        w  = Af / dot(rq,rq);

        M = M + w * (rq(:)*rq(:).');
        du = u(q,1) - u(p,1);
        dv = u(q,2) - u(p,2);
        bu = bu + w * du * rq(:);
        bv = bv + w * dv * rq(:);
    end

    trM = trace(M); if trM<=0, trM = 1; end
    M = M + (1e-12*trM)*eye(2);
    gradU(p,:) = (M\bu).';
    gradV(p,:) = (M\bv).';
end
end