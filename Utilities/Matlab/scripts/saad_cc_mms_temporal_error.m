%!/usr/bin/matlab
%Vanella - from McDermotts saad_mms_temporal_error.m
%02-20-2018
%saad_cc_mms_temporal_error.m

close all
clear all

plot_style

% Problem 1 parameters
EXPLICIT=1;
IMPLICIT=2;
r = 0.5;
nx = 512;
L = 2;
dx = L/nx;
x = -L/2+dx/2:dx:L/2-dx/2;
rho__0 = 5;
rho__1 = .5;
f = .5*(1+sin(2*pi*x/L));
rho = 1./((1-f)./rho__0 + f./rho__1);
% plot(x,rho)
% return
% T = 0.00390625; % total simulation time (CFL=1)
% dt = [T .5*T .25*T .125*T .0625*T];

% Output files

datadir = '../../Verification/Complex_Geometry/';
plotdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

for exim=EXPLICIT:IMPLICIT

   if exim == EXPLICIT
      filename = {'saad_CC_explicit_512_cfl_p25_mms.csv', ...
                  'saad_CC_explicit_512_cfl_p125_mms.csv', ...
                  'saad_CC_explicit_512_cfl_p0625_mms.csv'};
   else
      filename = {'saad_CC_implicit_512_cfl_p25_mms.csv', ...
                  'saad_CC_implicit_512_cfl_p125_mms.csv', ...
                  'saad_CC_implicit_512_cfl_p0625_mms.csv'};
   end

   skip_case = 0;

   for n=1:length(filename)
       if ~exist([datadir,filename{n}])
           display(['Error: File ' [datadir,filename{n}] ' does not exist. Skipping case.'])
           skip_case = 1;
       end
   end

   if skip_case
       return
   end

   % Gather FDS results

   n=0;
   M1 = importdata([datadir,filename{n+1}],',',2);
   M2 = importdata([datadir,filename{n+2}],',',2);
   M3 = importdata([datadir,filename{n+3}],',',2);

   ii = [nx+1:2*nx];
   rho_1 = M1.data(ii,1);
   rho_2 = M2.data(ii,1);
   rho_3 = M3.data(ii,1);
   Z_1 = M1.data(ii,2);
   Z_2 = M2.data(ii,2);
   Z_3 = M3.data(ii,2);

   p_rho = log( abs(rho_3-rho_2)./abs(rho_2-rho_1) )./log(r);
   p_Z = log( abs(Z_3-Z_2)./abs(Z_2-Z_1) )./log(r);

   % disp('Saad CC temporal order')
   % disp(' ')
   % disp(['L1 p rho = ',num2str( norm(p_rho,1)/nx )])
   % disp(['L2 p rho = ',num2str( norm(p_rho,2)/sqrt(nx) )])
   % disp(['Linf p rho = ',num2str( norm(p_rho,-inf) )])
   % disp(' ')
   % disp(['L1 p Z = ',num2str( norm(p_Z,1)/nx )])
   % disp(['L2 p Z = ',num2str( norm(p_Z,2)/sqrt(nx) )])
   % disp(['Linf p Z = ',num2str( norm(p_Z,-inf) )])
   % disp(' ')

   % flag errors

   L2_rho = norm(p_rho,2)/sqrt(nx);
   if L2_rho<1.99
       if exim == EXPLICIT
          disp(['Matlab Warning: L2_rho = ',num2str(L2_rho),' in Saad CC Explicit MMS'])
       else
          disp(['Matlab Warning: L2_rho = ',num2str(L2_rho),' in Saad CC Implicit MMS'])
       end
   end

   L2_Z = norm(p_Z,2)/sqrt(nx);
   if L2_Z<1.99
       if exim == EXPLICIT
          disp(['Matlab Warning: L2_Z = ',num2str(L2_Z),' in Saad CC Explicit MMS'])
       else
          disp(['Matlab Warning: L2_Z = ',num2str(L2_Z),' in Saad CC Implicit MMS'])
       end
   end

   % write the l2 norm to latex

   if exim == EXPLICIT
      fid = fopen([plotdir,'saad_CC_explicit_l2_norm.tex'],'wt');
   else
      fid = fopen([plotdir,'saad_CC_implicit_l2_norm.tex'],'wt');
   end
   fprintf(fid,'%s\n',num2str(L2_rho));
   fclose(fid);

   % Tony Saad's way...
   p1_rho_saad = log( norm(rho_3-rho_2,1)./norm(rho_2-rho_1,1) )./log(r);
   p2_rho_saad = log( norm(rho_3-rho_2,2)./norm(rho_2-rho_1,2) )./log(r);
   pinf_rho_saad = log( norm(rho_3-rho_2,inf)./norm(rho_2-rho_1,inf) )./log(r);
   % disp(['L1 p rho Saad = ',num2str( p1_rho_saad )])
   % disp(['L2 p rho Saad = ',num2str( p2_rho_saad )])
   % disp(['Linf p rho Saad = ',num2str( pinf_rho_saad )])
   % disp(' ')
   p1_Z_saad = log( norm(Z_3-Z_2,1)./norm(Z_2-Z_1,1) )./log(r);
   p2_Z_saad = log( norm(Z_3-Z_2,2)./norm(Z_2-Z_1,2) )./log(r);
   pinf_Z_saad = log( norm(Z_3-Z_2,inf)./norm(Z_2-Z_1,inf) )./log(r);
   % disp(['L1 p Z Saad = ',num2str( p1_Z_saad )])
   % disp(['L2 p Z Saad = ',num2str( p2_Z_saad )])
   % disp(['Linf p Z Saad = ',num2str( pinf_Z_saad )])
   % disp(' ')

   figure
   set(gca,'Units',Plot_Units)
   set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

   plot(x,p_rho)
   xlabel('{\it x} (m)','FontName',Font_Name,'FontSize',Label_Font_Size)
   ylabel('{\it p} order density','FontName',Font_Name,'FontSize',Label_Font_Size)
   if exim == EXPLICIT
      Git_Filename = [datadir,'saad_CC_explicit_512_cfl_p0625_git.txt'];
   else
      Git_Filename = [datadir,'saad_CC_implicit_512_cfl_p0625_git.txt'];
   end
   addverstr(gca,Git_Filename,'linear')

   set(gca,'FontName',Font_Name)
   set(gca,'FontSize',Label_Font_Size)
   set(gcf,'Visible',Figure_Visibility);
   set(gcf,'Units',Paper_Units);
   set(gcf,'PaperUnits',Paper_Units);
   set(gcf,'PaperSize',[Paper_Width Paper_Height]);
   set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
   if exim == EXPLICIT
      print(gcf,'-dpdf',[plotdir,'saad_CC_explicit_temporal_order_rho'])
   else
      print(gcf,'-dpdf',[plotdir,'saad_CC_implicit_temporal_order_rho'])
   end

   figure
   set(gca,'Units',Plot_Units)
   set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

   plot(x,p_Z)
   xlabel('{\it x} (m)','FontName',Font_Name,'FontSize',Label_Font_Size)
   ylabel('{\it p} order mixture fraction','FontName',Font_Name,'FontSize',Label_Font_Size)
   addverstr(gca,Git_Filename,'linear')

   set(gca,'FontName',Font_Name)
   set(gca,'FontSize',Label_Font_Size)
   set(gcf,'Visible',Figure_Visibility);
   set(gcf,'Units',Paper_Units);
   set(gcf,'PaperUnits',Paper_Units);
   set(gcf,'PaperSize',[Paper_Width Paper_Height]);
   set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
   if exim == EXPLICIT
      print(gcf,'-dpdf',[plotdir,'saad_CC_explicit_temporal_order_Z'])
   else
      print(gcf,'-dpdf',[plotdir,'saad_CC_implicit_temporal_order_Z'])
   end

   figure
   set(gca,'Units',Plot_Units)
   set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

   plot(x,rho,'r--'); hold on
   plot(x,rho_3,'r-')
   xlabel('{\it x} (m)','FontName',Font_Name,'FontSize',Label_Font_Size)
   ylabel('density (kg/m^3)','FontName',Font_Name,'FontSize',Label_Font_Size)
   lh=legend('initial field','final field','location','northeast');
   set(lh,'FontName',Font_Name,'FontSize',Label_Font_Size)
   addverstr(gca,Git_Filename,'linear')

   set(gca,'FontName',Font_Name)
   set(gca,'FontSize',Label_Font_Size)
   set(gcf,'Visible',Figure_Visibility);
   set(gcf,'Units',Paper_Units);
   set(gcf,'PaperUnits',Paper_Units);
   set(gcf,'PaperSize',[Paper_Width Paper_Height]);
   set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
   if exim == EXPLICIT
      print(gcf,'-dpdf',[plotdir,'saad_CC_explicit_rho'])
   else
      print(gcf,'-dpdf',[plotdir,'saad_CC_implicit_rho'])
   end

   figure
   set(gca,'Units',Plot_Units)
   set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

   plot(x,f,'b--'); hold on
   plot(x,Z_3,'b-')
   xlabel('{\it x} (m)','FontName',Font_Name,'FontSize',Label_Font_Size)
   ylabel('mixture fraction, {\it Z}','FontName',Font_Name,'FontSize',Label_Font_Size)
   lh=legend('initial field','final field','location','northwest');
   set(lh,'FontName',Font_Name,'FontSize',Label_Font_Size)
   addverstr(gca,Git_Filename,'linear')

   set(gca,'FontName',Font_Name)
   set(gca,'FontSize',Label_Font_Size)
   set(gcf,'Visible',Figure_Visibility);
   set(gcf,'Units',Paper_Units);
   set(gcf,'PaperUnits',Paper_Units);
   set(gcf,'PaperSize',[Paper_Width Paper_Height]);
   set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
   if exim == EXPLICIT
      print(gcf,'-dpdf',[plotdir,'saad_CC_explicit_Z'])
   else
      print(gcf,'-dpdf',[plotdir,'saad_CC_implicit_Z'])
   end

   figure
   set(gca,'Units',Plot_Units)
   set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

   plot(x,rho_3-rho_2,'b-'); hold on
   plot(x,rho_2-rho_1,'r-')
   xlabel('{\it x} (m)','FontName',Font_Name,'FontSize',Label_Font_Size)
   ylabel('density difference (kg/m^3) ','FontName',Font_Name,'FontSize',Label_Font_Size)
   lh=legend('\rho_3-\rho_2','\rho_2-\rho_1','location','northeast');
   set(lh,'FontName',Font_Name,'FontSize',Label_Font_Size)
   addverstr(gca,Git_Filename,'linear')

   set(gca,'FontName',Font_Name)
   set(gca,'FontSize',Label_Font_Size)
   set(gcf,'Visible',Figure_Visibility)
   set(gcf,'Units',Paper_Units)
   set(gcf,'PaperUnits',Paper_Units);
   set(gcf,'PaperSize',[Paper_Width Paper_Height]);
   set(gcf,'Position',[0 0 Paper_Width Paper_Height])
   if exim == EXPLICIT
      print(gcf,'-dpdf',[plotdir,'saad_CC_explicit_rho_diff'])
   else
      print(gcf,'-dpdf',[plotdir,'saad_CC_implicit_rho_diff'])
   end

end
