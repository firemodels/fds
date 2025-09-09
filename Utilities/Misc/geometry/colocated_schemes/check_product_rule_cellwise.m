function check_product_rule_cellwise(u_cell, H_cell, cells, faces)
Nc = numel(cells); Nf = numel(faces);
V  = reshape([cells.V],[],1);

% Green–Gauss grad H and Hf
sum_flux = zeros(Nc,2);
Hf = zeros(Nf,1);
for f = 1:Nf
    P  = faces(f).owner; 
    N  = faces(f).neigh; 
    nf = faces(f).nf(:); 
    Af = faces(f).Af;

    if N > 0
        Hf(f) = 0.5*(H_cell(P) + H_cell(N));
    else
        Hf(f) = H_cell(P);
    end
    contrib = (Hf(f) * nf.').*Af;  % 1x2
    sum_flux(P,:) = sum_flux(P,:) + contrib;
    if N > 0
        sum_flux(N,:) = sum_flux(N,:) - contrib;
    end
end
gradH = sum_flux ./ V;  % Nc x 2

% RHS and divU using cell-interpolated face-normal speed
rhs_cell = zeros(Nc,1);
divU     = zeros(Nc,1);
for p = 1:Nc
    S_rhs = 0; S_div = 0;
    fids = cells(p).faces;
    for k = 1:numel(fids)
        f  = fids(k);
        nf = faces(f).nf(:); 
        Af = faces(f).Af;
        s  = +1; if faces(f).neigh==p, s = -1; end

        P  = faces(f).owner; 
        N  = faces(f).neigh;
        if N > 0
            uf = 0.5*(u_cell(P,:) + u_cell(N,:));
        else
            uf = u_cell(P,:);
        end
        ufn = uf * nf;  % scalar

        S_rhs = S_rhs + s * Hf(f) * ufn * Af;
        S_div = S_div + s * ufn * Af;
    end
    rhs_cell(p) = S_rhs;
    divU(p)     = S_div / V(p);
end

lhs_cell = V .* ( sum(u_cell .* gradH, 2) + H_cell(:).*divU );
res_cell = rhs_cell - lhs_cell;

fprintf('Cell product rule: L2=%.3e  Linf=%.3e  sum(RHS)=%.3e  sum(LHS)=%.3e  diff=%.3e\n', ...
        norm(res_cell,2), norm(res_cell,inf), sum(rhs_cell), sum(lhs_cell), sum(rhs_cell)-sum(lhs_cell));
end
