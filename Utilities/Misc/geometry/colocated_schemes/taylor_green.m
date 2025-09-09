function [params,hist,cells,faces]=taylor_green(Nx_in,Ny_in,nu,dstep_plot_in)

% Nx; Ny  % cells in x,y
% nu       % kinematic viscosity
% Rest of parameters as in FDS Verification Guide Taylor-Green vortex
% problem (ns2d).
% Taylor green Vortex Problem, colocated grid, non-conservative Low Mach
% approximation:
% --------------------------------------------------------------
global compact_laplacian ux0 uy0 Nx Ny Xc Yc dstep_plot do_ssprk2 ...
       H_projection
%close all
%clc

compact_laplacian=1;
do_ssprk2 =1;
H_projection=1;

%% ------------------------ FDS params ------------------------
Nx      = Nx_in; Ny = Ny_in; % Number of cells inx and y.
Lx      = 2*pi; Ly = 2*pi;   % domain size
rho0    = 1.0;               % density
mu0     = rho0*nu;           % dynamic viscosity
U0      = 2.0;               % TG amplitude
kx      = 1;  ky   = 1;      % TG wavenumbers
ux0     = 1.;  uy0  = 1.;    % TG constant advection velocities in x,y.
gvec    = [0, 0];            % No gravity.
CFL     = 0.25;              % CFL constraint.
VN      = 0.5;               % Diffusion constraint.
t_begin = 0;                 % Time domain parameters..
t_end=2.*pi; 
dt = 0.001;
variable_time_step = 1;
dstep_plot=dstep_plot_in;

%% ------------------------ Build mesh -------------------------
% cells, faces, Uface (Nf×1), mu (Nc×1), rho (Nc×1), Lx, Ly
[cells, faces] = build_square_mesh_periodic_2d(Nx, Ny, Lx, Ly);
Nc = numel(cells);  Nf = numel(faces);
% Convenience vectors
xc = zeros(Nc,2);
V  = zeros(Nc,1);
for p = 1:Nc
    xc(p,:) = cells(p).xc;
    V(p)    = cells(p).V;
end

%% ------------------------ Fields init ------------------------
% Taylor–Green vortex (2D, periodic box)
u  = zeros(Nc,2); x  = xc(:,1);  y = xc(:,2);

% Coordinates 2d:
% Our linear index is id(i,j) = (j-1)*Nx + i, so reshape(...,[Nx,Ny])'.
Xc = reshape(xc(:,1), [Nx,Ny])';
Yc = reshape(xc(:,2), [Nx,Ny])';

rho = rho0 * ones(Nc,1);
mu  = mu0  * ones(Nc,1); lambda = -2/3 * mu;   % Stokes' hypothesis.

% Initial Fields:
% Face-normal velocities:
Uface = zeros(Nf,1);
% Do faces first then reconstruct_u_from_faces:
for f = 1:Nf
    P = faces(f).owner; N = faces(f).neigh; nf = faces(f).nf(:);
    xf= faces(f).cf;
    ubar(1) =  -U0 * cos(kx*xf(1)) .* sin(ky*xf(2)) + ux0;
    ubar(2) =   U0 * sin(kx*xf(1)) .* cos(ky*xf(2)) + uy0;
    Uface(f) = ubar * nf;   % dot product
end

% Cell centered velocities:
% u = -U0 cos(kx x) sin(ky y),  v = U0 sin(kx x) cos(ky y)
u(:,1) = -U0 * cos(kx*x) .* sin(ky*y) + ux0;
u(:,2) =  U0 * sin(kx*x) .* cos(ky*y) + uy0;

% Cell centered pressure (not used in P0 fractional step):
p = -(rho0*U0^2/4) * (cos(2*kx*x) + cos(2*ky*y));

%% Try SSPRK2 :
params = struct;
params.Lx = Lx;
params.Ly = Ly;
params.rho0 = rho0;
params.gvec = gvec;          % TG: no body force
params.CFL  = CFL;
params.VN   = VN;
params.divU_th = zeros(numel(cells),1);      % incompressible
params.t_begin = 0.0;
params.t_end   = t_end;
params.dt      = dt;
params.variable_time_step = variable_time_step;
params.max_steps = 1e5;
params.add_meanK = true;
params.poisson_opts = struct('fix','pin','pin_idx',1,'use_ilu',true);
params.tg_params=struct('U0',U0,'kx',kx,'ky',ky,'nu',nu,'ux0',ux0,'uy0',uy0);

% Cell at Nx/2, Ny/2:
% Choose the logical center
if mod(Nx,2)==0, i_mid = Nx/2; else, i_mid = (Nx+1)/2; end
if mod(Ny,2)==0, j_mid = Ny/2; else, j_mid = (Ny+1)/2; end
params.p_mid=find_cell_by_center_centered(cells,Nx,Ny,Lx,Ly,i_mid,j_mid);
fprintf('Center cell (i=%d,j=%d)->index p= %d\n',i_mid,j_mid,params.p_mid);

figure
[u,Uface,p,H_cell,t_end,hist]= run_ssprk2(u, Uface, p, rho, mu, lambda, ...
                          cells, faces, Lx, Ly, params);



return
end

% Some code snippets:

% Reconstruct cell centered u from face normal velocities (LS):
% u = reconstruct_u_from_faces_w(Uface, cells, faces);

% % Rescale velocities to match mean analytical KE:
% Vol = sum(V);
% rho_mean = sum(V .* rho(:)) / Vol;
% K_tgt = rho_mean * (U0^2/4)*Vol + (ux0^2+uy0^2)/2*Vol;
% [u_scaled, Uface_scaled, alpha, K_curr, K_tgt] = ...
%     rescale_to_target_ke(u, Uface, rho, cells, K_tgt);
% u = u_scaled;
% Uface = Uface_scaled;

% face normal velocities as averages of cell centered velocities:
% for f = 1:Nf
%     P = faces(f).owner; N = faces(f).neigh; nf = faces(f).nf(:);
%     if N > 0
%         ubar = 0.5*(u(P,:)+u(N,:));
%     else
%         ubar = zeros(1,2);
%     end
%     Uface(f) = ubar * nf;
% end