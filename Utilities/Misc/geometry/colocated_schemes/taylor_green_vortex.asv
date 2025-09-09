% taylor_green_vortex:
%
% Routine that calls Taylor-Green time integration for different number 
% of cells in Nx and Ny, as well as nu, and plots conservation, error norms 
% of colocated SSPRK2 solver.
% -------------------------------------------------------------------------
close all
clear all
clc

Nx=[8 16 32 64]; Ny=Nx;
nu=[0. 0.1];
dstep_plot = 100;

nNu = numel(nu); nNx = numel(Nx);
ErrM  = nan(nNu,nNx);          % L2 momentum abs error at T_end
RMSu  = nan(nNu,nNx);          % RMS u(mid) error
Nxv   = Nx(:)';                % row
for inu=1:length(nu)
    for ix=1:length(Nx)

        [params,hist,cells,faces]=taylor_green(Nx(ix),Ny(ix),nu(inu),dstep_plot);

        % Plot Mass (divergence), momentum and Kinetic Energy vs time:
        plot_mass_mom_ke(hist, cells);
        
        % Plot L2 errors of momentum and KE vs time:
        plot_L2_errors(hist);
        
        % --- grab last entry ---
        idx   = numel(hist.t);
        tend  = hist.t(idx);
        nsteps= hist.step(idx);
        
        L2_Mom_err   = hist.L2_m_abs(idx);
        L2_KE_err    = hist.L2_k_abs(idx);
        
        % --- print nicely ---
        fprintf('  nu, Nx, Ny =  %.6g, %s %.6g, %s %.6g, Final steps, time t = %5d, %s %.6g\n',nu(inu),' ',Nx(ix),' ',Ny(ix),nsteps,' ',tend);
        fprintf('  L2 momentum abs error   : %.3e\n', L2_Mom_err);
        fprintf('  L2 kinetic energy abs   : %.3e\n', L2_KE_err);
        
        % Now do L2 (RMS error) on the time series for velocities at p_mid:
        U0 = params.tg_params.U0;
        kx = params.tg_params.kx;
        ky = params.tg_params.ky;
        ux0= params.tg_params.ux0;
        uy0= params.tg_params.uy0;
        xp=cells(params.p_mid).xc(1); yp=cells(params.p_mid).xc(2);
        for i=1:idx
           t=hist.t(i);
           [up, vp] = taylor_green_uv_2d_point(xp,yp,t,U0,kx,ky,nu(inu),ux0,uy0);
           hist.u_ex_mid(i,1:2) = [up vp];
        end
        
        RMS_u_mid = sqrt(1/idx*sum((hist.u_num_mid(:,1)-hist.u_ex_mid(:,1)).^2));
        RMS_v_mid = sqrt(1/idx*sum((hist.u_num_mid(:,2)-hist.u_ex_mid(:,2)).^2));
        
        fprintf('  RMS u error mid   : %.3e\n', RMS_u_mid);
        fprintf('  RMS v error mid   : %.3e\n', RMS_v_mid);

        TG_Stats{inu,ix}.params = params;
        TG_Stats{inu,ix}.hist = hist;
        TG_Stats{inu,ix}.nsteps = nsteps;
        TG_Stats{inu,ix}.L2_Mom_err = L2_Mom_err;
        TG_Stats{inu,ix}.L2_KE_err  = L2_KE_err;
        TG_Stats{inu,ix}.RMS_u_mid  = RMS_u_mid;
        TG_Stats{inu,ix}.RMS_v_mid  = RMS_v_mid;

        % --- Plot u(midpoint) vs time: numerical vs exact ---
        t     = hist.t(:);
        u_num = hist.u_num_mid(:,1);   % numerical u at midpoint
        u_ex  = hist.u_ex_mid(:,1);    % exact u at midpoint
        
        figure('Name','Midpoint u vs time','Color','w');
        plot(t, u_num, '-o', 'LineWidth',1.8, 'MarkerSize',5, 'DisplayName','u_{num}(mid)');
        hold on
        plot(t, u_ex,  '--', 'LineWidth',2.0, 'DisplayName','u_{exact}(mid)');
        grid on
        xlabel('t','FontSize',14);
        ylabel('u at midpoint','FontSize',14);
        title(['Midpoint u(t): numerical vs exact ' num2str(Nx(ix)) '^2 cells'],'FontSize',16);
        legend('Location','best');
        set(gca,'FontSize',12);

    end

    % Plot L2 momentum error norm at T_end, and RMS error of u-velocity for
    % mid - point:
    h = TG_Stats{1,1}.params.Lx ./ Nxv;
    for ix = 1:nNx
        if ~isempty(TG_Stats{inu,ix})
            ErrM(inu,ix) = TG_Stats{inu,ix}.L2_Mom_err;
            RMSu(inu,ix) = TG_Stats{inu,ix}.RMS_u_mid;
        end
    end
    
    % ==== Plot L2 momentum error vs grid ====
    figure
    tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
    
    % L2 err momentum vs h = L_x / N_x
    nexttile;
    hold on
    plot(h, ErrM(inu,:), '-s','LineWidth',1.8,'MarkerSize',6, ...
            'DisplayName', sprintf('\\nu=%.3g', nu(inu)));
    plot(h,3*h,'--k','DisplayName', sprintf('O(\\Delta x)'))
    plot(h,0.8*h.^2,'-k','DisplayName', sprintf('O(\\Delta x^2)'))
    set(gca,'YScale','log','XScale','log','FontSize',14); grid on
    xlabel('\Delta x','FontSize',16); ylabel('||M||_2 at T_{end}','FontSize',16);
    title(['L2 momentum error vs h, nu=' num2str(nu(inu))],'FontSize',18);
    legend('Location','best','FontSize',16);
    axis([0.08 1. 10^-3 3])
    box on

    % ==== Plot RMS u(mid) error vs grid ====
    figure
    tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
        
    % RMS error mid point u-velocity vs h
    nexttile;
    hold on
    plot(h, RMSu(inu,:), '-s','LineWidth',1.8,'MarkerSize',6, ...
            'DisplayName', sprintf('Numerical \\nu=%.3g', nu(inu)));
    plot(h,2*h,'--k','DisplayName', sprintf('O(\\Delta x)'))
    plot(h,.3*h.^2,'-k','DisplayName', sprintf('O(\\Delta x^2)'))
    set(gca,'YScale','log','XScale','log','FontSize',14); grid on
    xlabel('\Delta x','FontSize',16); ylabel('U-velocity RMS error','FontSize',16);
    title(['U-velocity RMS error vs h, nu=' num2str(nu(inu))],'FontSize',18);
    legend('Location','best','FontSize',16);
    axis([0.08 1 10^-3 2])
    box on

end

return