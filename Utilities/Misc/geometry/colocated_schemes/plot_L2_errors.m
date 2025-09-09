function plot_L2_errors(hist)
%PLOT_L2_ERRORS  Plot L2 absolute errors for momentum and kinetic energy.
% Expects fields: hist.t, hist.L2_m_abs, hist.L2_k_abs (vectors)

t   = hist.t(:);
eMm = hist.L2_m_abs(:);
eKk = hist.L2_k_abs(:);

% Guard against zeros for log scale
eMm_plot = max(eMm, eps);
eKk_plot = max(eKk, eps);

figure('Name','L2 Errors','Color','w');

subplot(1,2,1)
semilogy(t, eMm_plot, '-o','LineWidth',1.8,'MarkerSize',5); grid on
xlabel('t','FontSize',16); ylabel('L2(M)','FontSize',16)
title('L2 momentum error','FontSize',18)
set(gca,'FontSize',14)

subplot(1,2,2)
semilogy(t, eKk_plot, '-s','LineWidth',1.8,'MarkerSize',5); grid on
xlabel('t','FontSize',16); ylabel('L2(KE)','FontSize',16)
title('L2 kinetic energy error','FontSize',18)
set(gca,'FontSize',14)

end
