function divU = cell_divergence_from_faces(Uface, cells, faces)
% Returns (∇·u)_p from updated face-normal velocities
Nc = numel(cells);
divU = zeros(Nc,1);
for p = 1:Nc
    S = 0.0;
    for k = 1:numel(cells(p).faces)
        f  = cells(p).faces(k);
        Af = faces(f).Af;
        s = +1; if faces(f).neigh==p, s = -1; end
        S = S + s * Uface(f) * Af;
    end
    divU(p) = S / cells(p).V;
end
end
