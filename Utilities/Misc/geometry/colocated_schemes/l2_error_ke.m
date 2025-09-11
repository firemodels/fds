function [L2_abs, L2_rel, K_num_tot, K_ex_tot, K_ex_formula] = l2_error_ke(u_num, rho, cells, t, tg)

global ux0 uy0

Nc = numel(cells);
xc = zeros(Nc,1); yc = zeros(Nc,1); V = zeros(Nc,1);
for p = 1:Nc, xc(p)=cells(p).xc(1); yc(p)=cells(p).xc(2); V(p)=cells(p).V; end
u_ex = taylor_green_uv_2d(xc, yc, t, tg.U0, tg.kx, tg.ky, tg.nu);
k_num = 0.5*rho(:).*sum(u_num.^2,2);
k_ex  = 0.5*rho(:).*sum(u_ex.^2,2);
num_sq = V.*(k_num-k_ex).^2;
den_sq = V.*(k_ex).^2;
L2_abs = sqrt(sum(num_sq)/sum(V));
L2_rel = sqrt(sum(num_sq)/max(1e-300,sum(den_sq)));
K_num_tot = sum(V.*k_num);
K_ex_tot  = sum(V.*k_ex);
% Correct analytic total KE (mean k = ρ U0^2 / 4 e^{-2 ν k^2 t})
k2 = tg.kx*tg.kx + tg.ky*tg.ky;
Vol = sum(V);
rho_bar = sum(V.*rho(:))/Vol;
K_ex_formula = rho_bar*(tg.U0^2/4)*exp(-2*tg.nu*k2*t)*Vol+(ux0^2+uy0^2)/2*Vol;
end
