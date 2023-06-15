% Mueller
% 06-15-2023
% tree_shape.m

% This script compares the actual mass generated with various tree crown
% geometries to the expected value (assuming a bulk density of 5 kg/m3)

close all
clear all

out_dir = '../../Verification/WUI/';
plot_dir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

% import mass totals
CHID = [out_dir,'tree_shapes'];
if ~exist([out_dir,CHID,'_devc.csv'],'file') 
    disp(['Error: File ',CHID,'_devc.csv does not exist. Skipping case.']);
    return;
end
DEVC = importdata([out_dir,[CHID,'_devc.csv']],',', 2);

random_part_mass=DEVC.data(end,2:6);
one_part_mass=DEVC.data(end,7:end);


% compute ideal input volumes
box_vol = 0.6*0.6*0.4;
cone_vol = 0.6*0.3^2*pi/3;
cyl_vol = 0.6*0.2^2*pi;
cone_shell_vol = 0.6*(0.3^2-0.2^2)*pi/3;
cyl_shell_vol = 0.6*(0.2^2-0.15^2)*pi;

input_vol=[box_vol,cone_vol,cyl_vol,cone_shell_vol,cyl_shell_vol];

% relative error in input volume
% large tolerance because error in mass does not necessarily imply an issue
% as long as MASS_PER_VOLUME is preserved
rel_err=abs([random_part_mass one_part_mass]-5*[input_vol input_vol])./...
    (5*[input_vol input_vol]);
if max(rel_err) > 0.1
   display(['Matlab Warning: The mass in tree_shapes is out of tolerance.'])
end

% figure setup
figure
plot_style
colors=['k','r','b','g'];
clf
hold on
box on 

% plot
plot([0,.2],[0,1],'k--','displayname','ideal');
plot(input_vol,random_part_mass,'ro','displayname','1000 random particles');
plot(input_vol,one_part_mass,'bo','displayname','1 particle per cell');

% standardize figure
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Scat_Label_Font_Size)    
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

% label axes
axis([0 .15 0 .8])
xlabel('Input volume (m^3)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Tree crown mass (kg)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)     
% add legend
lh=legend();
set(lh,'FontSize',Key_Font_Size,'location','southeast')
% add Git revision
Git_Filename = [out_dir,CHID,'_git.txt'];
addverstr(gca,Git_Filename,'linear')
% save figure
print(gcf,Image_File_Type,[plot_dir,'tree_shapes'])

