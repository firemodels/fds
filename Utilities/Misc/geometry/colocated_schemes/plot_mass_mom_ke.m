function plot_mass_mom_ke(hist, cells)
%PLOT_MASS_MOM_KE  Plot mass residual, <momentum>, and <kinetic energy>.
% Inputs:
%   hist  : struct from run_ssprk2 (fields used: t, L2_mass, M_num, M_ex,
%           K_num_tot, K_ex_formula)
%   cells : mesh cells (used to compute domain volume)

Vol = sum([cells.V]);

% Averages
Mx_num = hist.M_num(:,1) / Vol;
My_num = hist.M_num(:,2) / Vol;
Mmag_num = sqrt(Mx_num.^2 + My_num.^2);

hasAnalyticM = isfield(hist,'M_ex') && all(size(hist.M_ex)==size(hist.M_num)) ...
               && all(~isnan(hist.M_ex(:)));
if hasAnalyticM
    Mx_ex = hist.M_ex(:,1) / Vol;  My_ex = hist.M_ex(:,2) / Vol;
    Mmag_ex = sqrt(Mx_ex.^2 + My_ex.^2);
else
    Mx_ex = []; My_ex = []; Mmag_ex = [];
end

KE_num_avg = hist.K_num_tot / Vol;
hasAnalyticK = isfield(hist,'K_ex_formula') && all(~isnan(hist.K_ex_formula(:)));
if hasAnalyticK
    KE_ex_avg = hist.K_ex_tot / Vol; %hist.K_ex_formula 
else
    KE_ex_avg = [];
end

% ---- Plots ----
figure('Name','Mass/Momentum/KE','Color','w','Position',[100 100 1200 400]);

fs_axis  = 12;   % axis/tick font size
fs_label = 14;   % x/y label size
fs_title = 16;   % title size
fs_leg   = 14;   % legend size

% (1) Mass conservation residual
subplot(1,3,1)
semilogy(hist.t, hist.L2_mass, '-o','LineWidth',1.5,'MarkerSize',6); grid on
xlabel('t','FontSize',fs_label,'FontWeight','bold');
ylabel('||mass||_2','FontSize',fs_label,'FontWeight','bold');
title('Mass conservation residual','FontSize',fs_title,'FontWeight','bold');
set(gca,'FontSize',fs_axis);

% (2) Average momentum
subplot(1,3,2)
plot(hist.t, Mmag_num, '-o', 'LineWidth',1.5,'MarkerSize',6,'DisplayName','|<M>| num'); hold on
plot(hist.t, Mx_num, '-s', 'LineWidth',1.5,'MarkerSize',6,'DisplayName','<M_x> num');
plot(hist.t, My_num, '-d', 'LineWidth',1.5,'MarkerSize',6,'DisplayName','<M_y> num');
if ~isempty(Mmag_ex)
    plot(hist.t, Mmag_ex, '--','LineWidth',1.5,'DisplayName','|<M>| ex');
    plot(hist.t, Mx_ex, '--','LineWidth',1.5,'DisplayName','<M_x> ex');
    plot(hist.t, My_ex, '--','LineWidth',1.5,'DisplayName','<M_y> ex');
end
hold off; grid on
xlabel('t','FontSize',fs_label,'FontWeight','bold');
ylabel('Momentum / Volume','FontSize',fs_label,'FontWeight','bold');
title('Domain-averaged momentum','FontSize',fs_title,'FontWeight','bold');
legend('Location','best','FontSize',fs_leg);
set(gca,'FontSize',fs_axis);

% (3) Average kinetic energy
subplot(1,3,3)
plot(hist.t, KE_num_avg, '-o','LineWidth',1.5,'MarkerSize',6,'DisplayName','<KE> num'); hold on
if ~isempty(KE_ex_avg)
    plot(hist.t, KE_ex_avg, '--','LineWidth',1.5,'DisplayName','<KE> ex');
end
hold off; grid on
xlabel('t','FontSize',fs_label,'FontWeight','bold');
ylabel('KE / Volume','FontSize',fs_label,'FontWeight','bold');
title('Domain-averaged kinetic energy','FontSize',fs_title,'FontWeight','bold');
legend('Location','best','FontSize',fs_leg);
set(gca,'FontSize',fs_axis);

sgtitle('Diagnostics: Volume-averaged quantities','FontSize',fs_title+2,'FontWeight','bold');

end