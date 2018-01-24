% McDermott & Hemley
% 6-26-12
% pressure_coefficient.m

close all
clear all

plot_style

res_dir = '../../../out/Wind_Engineering/FDS_Output_Files/';
exp_dir = '../../../exp/Wind_Engineering/';
plt_dir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Wind_Engineering/';

str_angle = {'180','270'}; % angles used
str_mesh = {'32','64'}; % cell numbers used
line_pos = {'-0p0009y','-0p0296y','-0p0437y','0p0256z','-0p0067y','-0p0368y','-0p0669y','0p0254z';... % 1st row: positions in filenames
    '-0.0009','-0.0296','-0.0437','0.0256','-0.0067','-0.0368','-0.0669','0.0254'}; % 2nd row: positions for plots
lines_per_angle = length(line_pos)/length(str_angle);

height = 0.0396; % building height in meters

for ind1 = 1:length(str_mesh) % cycles through meshes

    mesh = str_mesh{ind1};

    for ind2 = 1:length(str_angle) % cycles through angles

        col = 1;
        angle = str_angle{ind2};

        for ind3 = 1:lines_per_angle % cycles through lines

            line_file = line_pos{1,(ind3+4*(ind2-1))};
            line_plot = line_pos{2,(ind3+4*(ind2-1))};

            M = importdata([res_dir,'UWO_test7_case1_',angle,'_',mesh,'_line.csv'],',',2);
            N = importdata([exp_dir,'UWO_exp_',angle,'_',line_file,'.csv'],',',1);

            x_exp = N.data(:,1);
            cp_mean_exp = N.data(:,2);
            cp_rms_exp = N.data(:,3);

            if ind3 < lines_per_angle

                x = cell(3,1);
                cp = cell(3,1);

                for ind4 = 1:3 % cycles from windward wall to roof to leeward wall

                    x{ind4} = M.data(:,col);
                    cp{ind4} = M.data(:,col+1);

                    col = col+4;

                end

                x_1 = x{1};
                x_2 = height + x{2} - min(x{2});
                x_3 = max(x_2) + height - x{3};

                x_max_plot = ceil(max([x_3;x_exp])*100)/100;

                set(gca,'Units',Plot_Units)
                set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
                plot(x_1,cp{1},'bo-'); hold on
                plot(x_2,cp{2},'go-');
                plot(x_3,cp{3},'ro-');
                errorbar(x_exp,cp_mean_exp,cp_rms_exp,'ko');
                hold off
                set(gca,'FontName',Font_Name)
                set(gca,'FontSize',Label_Font_Size)
                axis([0,x_max_plot,-1.7,1.5])
                xlabel('Position along line (m)')
                ylabel('Mean C_p')
                Plot_Title=['{\it y} = ',line_plot,' m'];
                text(.01,1.2,Plot_Title,'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)

                % Add git version if file is available
                git_file = [res_dir,'UWO_test7_case1_',angle,'_',mesh,'_git.txt'];
                addverstr(gca,git_file,'linear')

                lh=legend('FDS (Windward Wall)','FDS (Roof)','FDS (Leeward Wall)', 'Exp','Location', 'NorthEast');
                set(lh,'FontSize',Key_Font_Size)
                set(lh,'Interpreter',Font_Interpreter)

                set(gcf,'Visible',Figure_Visibility);
                set(gcf,'Units',Paper_Units);
                set(gcf,'PaperUnits',Paper_Units);
                set(gcf,'PaperSize',[Paper_Width Paper_Height]);
                set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
                print('-dpdf',[plt_dir,'cp_mean_',angle,'_',mesh,'_',line_file])

            else
                % side wall
                x = M.data(:,col);
                x = x - min(x);
                cp = M.data(:,col+1);

                x_max_plot = ceil(max([x;x_exp])*100)/100;

                set(gca,'Units',Plot_Units)
                set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
                plot(x,cp,'bo-'); hold on
                errorbar(x_exp,cp_mean_exp,cp_rms_exp,'ko');
                hold off
                set(gca,'FontName',Font_Name)
                set(gca,'FontSize',Label_Font_Size)
                axis([0,x_max_plot,-1.7,1.5])
                xlabel('Position along line (m)')
                ylabel('Mean C_p')
                % title(['Mean C_p for an angle of ',angle,...
                %     ' on side wall at z = ',line_plot,' m'])

                Plot_Title=['Side wall, {\it z} = ',line_plot,' m'];
                text(.005,1.2,Plot_Title,'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)

                % Add git version if file is available
                git_file = [res_dir,'UWO_test7_case1_',angle,'_',mesh,'_git.txt'];
                addverstr(gca,git_file,'linear')

                lh=legend('FDS','Exp','Location', 'NorthEast');
                set(lh,'FontSize',Key_Font_Size)
                set(lh,'Interpreter',Font_Interpreter)

                set(gcf,'Visible',Figure_Visibility);
                set(gcf,'Units',Paper_Units);
                set(gcf,'PaperUnits',Paper_Units);
                set(gcf,'PaperSize',[Paper_Width Paper_Height]);
                set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
                print('-dpdf',[plt_dir,'cp_mean_',angle,'_',mesh,'_',line_file])

            end

            clear x x_1 x_2 x_3 cp M N x_exp cp_mean_exp cp_rms_exp

        end

    end

end

