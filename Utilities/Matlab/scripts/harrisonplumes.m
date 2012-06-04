%%harrisonplumes.m
%Taylor Myers
%Harrison_Spill_Plumes scatterplot maker
clear all; close all;

addpath '../../Validation/Harrison_Spill_Plumes/Experimental_Data'
addpath '../../Validation/Harrison_Spill_Plumes/FDS_Output_Files'

%Read in Experimental data
HSPstruct = importdata('Harrison_Spill_Plumes.csv'); %read in csv file
A=HSPstruct.data; %numerical data
At=HSPstruct.textdata; %text data
C=At(:,1:2);
S = zeros(length(A),2);
S(1:length(A),1) = A(1:length(A),4);

%Read in FDS data
%There was no neat way to do this. Forgive me for what I am about to do.
%SE4
Bstruct = importdata('SE4_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(7,2)  =  B(x,3);
S(10,2) =  B(x,3);
S(43,2) =  B(x,4);
S(46,2) =  B(x,4);
S(58,2) =  B(x,5);
S(61,2) =  B(x,5);
S(91,2) =  B(x,6);
S(94,2) =  B(x,6);
S(115,2)=  B(x,7);
S(118,2)=  B(x,7);
%SE5
Bstruct = importdata('SE5_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(8,2)  =  B(x,3);
S(11,2) =  B(x,3);
S(44,2) =  B(x,4);
S(47,2) =  B(x,4);
S(59,2) =  B(x,5);
S(62,2) =  B(x,5);
S(92,2) =  B(x,6);
S(95,2) =  B(x,6);
S(116,2)=  B(x,7);
S(119,2)=  B(x,7);
%SE6
Bstruct = importdata('SE6_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(9,2)  =  B(x,3);
S(12,2) =  B(x,3);
S(45,2) =  B(x,4);
S(60,2) =  B(x,5);
S(93,2) =  B(x,6);
S(117,2)=  B(x,7);
%SE7
Bstruct = importdata('SE7_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output;
[x,y]= size(B);
S(13,2) =  B(x,3);
S(16,2) =  B(x,3);
S(37,2) =  B(x,4);
S(40,2) =  B(x,4);
S(63,2) =  B(x,5);
S(66,2) =  B(x,5);
S(96,2) =  B(x,6);
S(99,2) =  B(x,6);
S(120,2)=  B(x,7);
S(123,2)=  B(x,7);
%SE8
Bstruct = importdata('SE8_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(14,2) =  B(x,3);
S(17,2) =  B(x,3);
S(38,2) =  B(x,4);
S(41,2) =  B(x,4);
S(64,2) =  B(x,5);
S(67,2) =  B(x,5);
S(97,2) =  B(x,6);
S(100,2)=  B(x,6);
S(121,2)=  B(x,7);
S(124,2)=  B(x,7);
%SE9
Bstruct = importdata('SE9_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(15,2) =  B(x,3);
S(18,2) =  B(x,3);
S(39,2) =  B(x,4);
S(42,2) =  B(x,4);
S(65,2) =  B(x,5);
S(98,2) =  B(x,6);
S(122,2)=  B(x,7);
%SE10
Bstruct = importdata('SE10_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(1,2)  =  B(x,3);
S(4,2)  =  B(x,3);
S(48,2) =  B(x,4);
S(51,2) =  B(x,4);
S(53,2) =  B(x,5);
S(56,2) =  B(x,5);
S(86,2) =  B(x,6);
S(89,2) =  B(x,6);
S(110,2)=  B(x,7);
S(113,2)=  B(x,7);
%SE11
Bstruct = importdata('SE11_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(2,2)  =  B(x,3);
S(5,2)  =  B(x,3);
S(49,2) =  B(x,4);
S(52,2) =  B(x,4);
S(54,2) =  B(x,5);
S(57,2) =  B(x,5);
S(87,2) =  B(x,6);
S(90,2) =  B(x,6);
S(111,2)=  B(x,7);
S(114,2)=  B(x,7);
%SE12
Bstruct = importdata('SE12_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(3,2)  =  B(x,3);
S(6,2)  =  B(x,3);
S(50,2) =  B(x,4);
S(55,2) =  B(x,5);
S(88,2) =  B(x,6);
S(112,2)=  B(x,7);
%SE13
Bstruct = importdata('SE13_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(22,2) =  B(x,3);
S(31,2) =  B(x,4);
S(71,2) =  B(x,5);
S(80,2) =  B(x,6);
S(104,2)=  B(x,7);
%SE14
Bstruct = importdata('SE14_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(23,2) =  B(x,3);
S(32,2) =  B(x,4);
S(72,2) =  B(x,5);
S(81,2) =  B(x,6);
S(105,2)=  B(x,7);
%SE15
Bstruct = importdata('SE15_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(24,2) =  B(x,3);
S(33,2) =  B(x,4);
S(73,2) =  B(x,5);
S(82,2) =  B(x,6);
S(106,2)=  B(x,7);
%SE16
Bstruct = importdata('SE16_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(19,2) =  B(x,2);
S(34,2) =  B(x,3);
S(68,2) =  B(x,4);
S(83,2) =  B(x,5);
S(101,2)=  B(x,6);
%SE17
Bstruct = importdata('SE17_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(20,2) =  B(x,2);
S(35,2) =  B(x,3);
S(69,2) =  B(x,4);
S(84,2) =  B(x,5);
S(102,2)=  B(x,6);
%SE18
Bstruct = importdata('SE18_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(21,2) =  B(x,2);
S(36,2) =  B(x,3);
S(70,2) =  B(x,4);
S(85,2) =  B(x,5);
S(103,2)=  B(x,6);
%SE19
Bstruct = importdata('SE19_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(25,2) =  B(x,3);
S(28,2) =  B(x,4);
S(74,2) =  B(x,5);
S(77,2) =  B(x,6);
S(107,2)=  B(x,7);
%SE20
Bstruct = importdata('SE20_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(26,2) =  B(x,3);
S(29,2) =  B(x,4);
S(75,2) =  B(x,5);
S(78,2) =  B(x,6);
S(108,2)=  B(x,7);
%SE21
Bstruct = importdata('SE21_devc.csv');%read in csv file
B=Bstruct.data; %numerical data from FDS output
[x,y]= size(B);
S(27,2) =  B(x,3);
S(30,2) =  B(x,4);
S(76,2) =  B(x,5);
S(79,2) =  B(x,6);
S(109,2)=  B(x,7);

%Format and Print Scatter Plot

paper_width  = 6.0; % inches
paper_height = 6.0; % inches
Plot_Min = 0; 
Plot_Max = .7; 
Sigma_E = .1;
Ind_Title = 'Measured Mass Flow kg/s';
Dep_Title = 'Predicted Mass Flow kg/s';
Scatter_Plot_Title = 'Harrison Spill Plumes Mass Flows';
Title_Position = [0.03 0.95];
Plot_Filename = 'FDS_Validation_Guide/FIGURES/ScatterPlots/Harrison_Spill_Plumes';

close all;
figure(1)
plot_style
hold on
plot(-1,-1, 'ko')
plot(-1,-1, 'r^')
for kk = 2:length(C)
    TF = strcmp(C(kk,2),'Balcony');
    if TF == 1
        plot(S(kk-1,1),S(kk-1,2),'ko');
    else
        plot(S(kk-1,1),S(kk-1,2),'r^');
    end
end
plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max],'k-')                    
plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max*(1+2*Sigma_E)],'k--') 
plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max*(1-2*Sigma_E)],'k--') 

% format the legend and axis labels
xlabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',14)
ylabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',14)
axis([Plot_Min Plot_Max Plot_Min Plot_Max])
legend('3D Balcony','3D Adhered','Location','SouthEast')


set(gca,'Units','inches')
set(gca,'FontName','Times')
set(gca,'FontSize',12)
set(gca,'YTick',get(gca,'XTick'))
set(gca,'Position',[1,1,4.5,4.5])
        
text(Plot_Min+Title_Position(1)*(Plot_Max-Plot_Min),Plot_Min+Title_Position(2)*(Plot_Max-Plot_Min),...
            Scatter_Plot_Title,'FontSize',14,'FontName','Times','Interpreter',Font_Interpreter)
  
if Sigma_E > 0.0
    text(Plot_Min+(Title_Position(1)+0.05)*(Plot_Max-Plot_Min),Plot_Min+(Title_Position(2)-0.05)*(Plot_Max-Plot_Min),...
                 ['$2 \, \sigma_E$=',num2str(2*Sigma_E,'%4.2f')],'FontSize',12,'FontName','Times','Interpreter',Font_Interpreter)
end
hold off
% print to pdf
        set(gcf,'Visible','on');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperSize',[paper_width paper_height]);
        set(gcf,'PaperPosition',[0 0 paper_width paper_height]);
        display(['Printing scatter plot ',num2str(j),'...'])
        print(gcf,'-dpdf',[pwd,'/../../Manuals/',Plot_Filename])


display('harrisonplumes completed successfully!')


