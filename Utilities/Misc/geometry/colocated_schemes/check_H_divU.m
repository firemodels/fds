% function check_H_divU(H_cell, Uface, cells, faces)
% Nc = numel(cells); V = reshape([cells.V],[],1);
% 
% % cell divergence from projected face speeds
% divU = zeros(Nc,1);
% for p=1:Nc
%   S=0; for k=1:numel(cells(p).faces)
%     f=cells(p).faces(k); s=+1; if faces(f).neigh==p, s=-1; end
%     S = S + s*Uface(f)*faces(f).Af;
%   end
%   divU(p) = S / V(p);
% end
% 
% % face-averaged H
% Nf = numel(faces); Hf = zeros(Nf,1);
% for f=1:Nf
%   P=faces(f).owner; N=faces(f).neigh;
%   Hf = (N > 0) * 0.5*(H_cell(P) + H_cell(N)) + (N == 0) * H_cell(P);
% end
% 
% LHS = sum( V .* (H_cell(:).*divU) );
% RHS = sum( Hf .* Uface .* [faces.Af].' );
% fprintf('(H,divU) face pairing: LHS=%.3e  RHS=%.3e  diff=%.3e\n', LHS, RHS, LHS-RHS);
% end

function check_H_divU(H_cell, Uface, cells, faces)
Nc = numel(cells); V = reshape([cells.V],[],1);

% cell divergence from faces (owner→neigh, s_pf = +1 for owner, -1 for neigh)
divU = zeros(Nc,1);
for p=1:Nc
    acc = 0.0;
    flist = cells(p).faces(:);
    for k=1:numel(flist)
        f = flist(k);
        s = +1; if faces(f).neigh==p, s = -1; end
        acc = acc + s * Uface(f) * faces(f).Af;
    end
    divU(p) = acc / V(p);
end

LHS = sum( V .* (H_cell(:) .* divU) );

% RHS = sum_f (H_owner - H_neigh) * U_f * A_f  (internal faces only)
owners = [faces.owner]';  neighs = [faces.neigh]';  Af = [faces.Af]';
mask = (neighs > 0);
dH = H_cell(owners(mask)) - H_cell(neighs(mask));
RHS = sum( dH .* Uface(mask) .* Af(mask) );

fprintf('(H,divU) pure pairing: LHS=%.6e  RHS=%.6e  diff=%.3e\n', LHS, RHS, LHS-RHS);
end
