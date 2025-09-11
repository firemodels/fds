function [L2_abs, L2_rel, M_num, M_ex] = l2_error_momentum(u_num, rho, cells, t, tg)
% tg has fields: U0, kx, ky, nu
Nc = numel(cells);
xc = zeros(Nc,1); yc = zeros(Nc,1); V = zeros(Nc,1);
for p = 1:Nc, xc(p)=cells(p).xc(1); yc(p)=cells(p).xc(2); V(p)=cells(p).V; end
u_ex = taylor_green_uv_2d(xc, yc, t, tg.U0, tg.kx, tg.ky, tg.nu);
m_num = rho(:).*u_num;  m_ex = rho(:).*u_ex;
num_sq = V.*sum((m_num-m_ex).^2,2);
den_sq = V.*sum(m_ex.^2,2);
L2_abs = sqrt(sum(num_sq)/sum(V));
L2_rel = sqrt(sum(num_sq)/max(1e-300,sum(den_sq)));
M_num = [sum(V.*m_num(:,1)), sum(V.*m_num(:,2))];
M_ex  = [sum(V.*m_ex(:,1)),  sum(V.*m_ex(:,2)) ];
end