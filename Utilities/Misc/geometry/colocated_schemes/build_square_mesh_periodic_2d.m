function [cells, faces] = build_square_mesh_periodic_2d(Nx, Ny, Lx, Ly)
    dx = Lx/Nx;  dy = Ly/Ny;
    Nc = Nx*Ny;
    xw = -Lx/2;
    ys = -Ly/2;
    
    % --- Cells ---
    cells(Nc) = struct('xc',[],'V',[],'faces',[]);
    id = @(i,j) (j-1)*Nx + i;             % 1-based linear index
    for j = 1:Ny
        yc = ys + (j-0.5)*dy;
        for i = 1:Nx
            xc_ij = xw + (i-0.5)*dx;
            p = id(i,j);
            cells(p).xc = [xc_ij, yc];
            cells(p).V  = dx*dy;
            cells(p).faces = [];          % will fill later
        end
    end

    % --- Faces (right and top only; periodic wraps add neighbors) ---
    faces = struct('owner',{},'neigh',{},'nf',{},'Af',{},'cf',{});
    fcount = 0;

    % 1) Vertical faces: between (i,j) and (i+1,j), wrap at i=Nx → 1
    for j = 1:Ny
        yc = ys + (j-0.5)*dy;
        for i = 1:Nx
            P = id(i,j);
            ir = i+1; if ir>Nx, ir = 1; end
            N = id(ir,j);

            cf_x = xw + i*dx; if cf_x>Lx, cf_x = cf_x - Lx; end
            cf = [cf_x, yc];
            nf = [1, 0];          % normal from P to N
            Af = dy;

            fcount = fcount + 1;
            faces(fcount).owner = P;
            faces(fcount).neigh = N;
            faces(fcount).nf    = nf;
            faces(fcount).Af    = Af;
            faces(fcount).cf    = cf;

            cells(P).faces(end+1) = fcount;
            cells(N).faces(end+1) = fcount;
        end
    end

    % 2) Horizontal faces: between (i,j) and (i,j+1), wrap at j=Ny → 1
    for j = 1:Ny
        jr = j+1; if jr>Ny, jr = 1; end
        y_face = ys + j*dy; if y_face>Ly, y_face = y_face - Ly; end
        for i = 1:Nx
            P = id(i,j);
            Np = id(i,jr);

            xc = xw + (i-0.5)*dx;
            cf = [xc, y_face];
            nf = [0, 1];          % normal from P to N
            Af = dx;

            fcount = fcount + 1;
            faces(fcount).owner = P;
            faces(fcount).neigh = Np;
            faces(fcount).nf    = nf;
            faces(fcount).Af    = Af;
            faces(fcount).cf    = cf;

            cells(P).faces(end+1) = fcount;
            cells(Np).faces(end+1) = fcount;
        end
    end

    % Build dpn:
    for f = 1:numel(faces)
    P = faces(f).owner; N = faces(f).neigh;
    nf = faces(f).nf(:);
    if N > 0
        dr_raw = cells(N).xc(:) - cells(P).xc(:);
        dr = periodic_delta_2d(dr_raw, Lx, Ly);
        faces(f).dPN = abs( dr' * nf );
    else
        faces(f).dPN = [];  % or a BC-specific value
    end
end

end

