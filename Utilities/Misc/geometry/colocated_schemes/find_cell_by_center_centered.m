function p = find_cell_by_center_centered(cells, Nx, Ny, Lx, Ly, i, j)
% Find the cell index p whose centroid is closest to logical indices (i,j)
% on a uniform grid spanning [-Lx/2, Lx/2] × [-Ly/2, Ly/2] with periodic wrap.

dx = Lx / Nx; 
dy = Ly / Ny;

% target cell-center coordinates for (i,j)
xt = -Lx/2 + (i - 0.5) * dx;
yt = -Ly/2 + (j - 0.5) * dy;

% collect centroids
Nc  = numel(cells);
Xc  = zeros(Nc,2);
for k = 1:Nc
    Xc(k,:) = cells(k).xc(:).';
end

% minimal-image periodic distances to target (works for centered domain)
dxv = Xc(:,1) - xt;  dxv = dxv - round(dxv / Lx) * Lx;
dyv = Xc(:,2) - yt;  dyv = dyv - round(dyv / Ly) * Ly;

[~, p] = min(dxv.^2 + dyv.^2);
end