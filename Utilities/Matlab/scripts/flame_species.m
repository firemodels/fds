%------------------------------------------------
% C Weinschenk
% 10-2011
% Combine outputs from Methane_flame_lumped and
% Methane_flame_primitive into 1 file
%------------------------------------------------

%-----------------------
% Write directory
%-----------------------

FDS_Output_Files = '../../Verification/Species/';

%----------------------
% Import files
%----------------------

if ~exist([FDS_Output_Files,'methane_flame_primitive_devc.csv'])
    display(['Error: File ',[FDS_Output_Files,'methane_flame_primitive_devc.csv'],' does not exist. Skipping case.'])
    return
end

if ~exist([FDS_Output_Files,'methane_flame_lumped_devc.csv'])
    display(['Error: File ',[FDS_Output_Files,'methane_flame_lumped_devc.csv'],' does not exist. Skipping case.'])
    return
end

if ~exist([FDS_Output_Files,'methane_flame_lumped_fuel_devc.csv'])
    display(['Error: File ',[FDS_Output_Files,'methane_flame_lumped_fuel_devc.csv'],' does not exist. Skipping case.'])
    return
end

if ~exist([FDS_Output_Files,'methane_flame_lumped_ox_devc.csv'])
    display(['Error: File ',[FDS_Output_Files,'methane_flame_lumped_ox_devc.csv'],' does not exist. Skipping case.'])
    return
end

if ~exist([FDS_Output_Files,'methane_flame_primitive_2_devc.csv'])
    display(['Error: File ',[FDS_Output_Files,'methane_flame_primitive_2_devc.csv'],' does not exist. Skipping case.'])
    return
end


primitive_struct=importdata([FDS_Output_Files,'methane_flame_primitive_devc.csv']);
primitive_data=primitive_struct.data;
primitive_text=primitive_struct.textdata;

lumped_struct=importdata([FDS_Output_Files,'methane_flame_lumped_devc.csv']);
lumped_data=lumped_struct.data;
lumped_text=lumped_struct.textdata;

header1=cat(2,primitive_text(1),lumped_text(1,2:4));
header2=cat(2,primitive_text(2,:),lumped_text(2,2:4));
data=cat(2,primitive_data,lumped_data(:,2:4));

primitive_struct_2=importdata([FDS_Output_Files,'methane_flame_primitive_2_devc.csv']);
primitive_data_2=primitive_struct_2.data;
primitive_text_2=primitive_struct_2.textdata;

lumped_struct_f=importdata([FDS_Output_Files,'methane_flame_lumped_fuel_devc.csv']);
lumped_data_f=lumped_struct_f.data;
lumped_text_f=lumped_struct_f.textdata;

lumped_struct_ox=importdata([FDS_Output_Files,'methane_flame_lumped_ox_devc.csv']);
lumped_data_ox=lumped_struct_ox.data;
lumped_text_ox=lumped_struct_ox.textdata;

header1b=cat(2,primitive_text_2(1),lumped_text_f(1,2:3),lumped_text_ox(1,2:3));
header2b=cat(2,primitive_text_2(2,:),lumped_text_f(2,2:3),lumped_text_ox(2,2:3));
datab=cat(2,primitive_data_2,lumped_data_f(:,2:3),lumped_data_ox(:,2:3));

%-----------------------
% Write new files
%-----------------------

fid = fopen([FDS_Output_Files,'methane_flame_lumpedprimitive.csv'],'wt','n');

fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n','s','kg','kg','kg','kg','kg','kg');
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n',header2{1,:});
for i=1:length(primitive_data)
    fprintf(fid,'%f,%f,%f,%f,%f,%f,%f\n',data(i,:));
end
fclose(fid);

fid = fopen([FDS_Output_Files,'methane_flame_multilumped.csv'],'wt','n');

fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n','s','kg','kg','kg','kg','kg','kg');
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n',header2b{1,:});
for i=1:length(primitive_data_2)
    fprintf(fid,'%f,%f,%f,%f,%f,%f,%f\n',datab(i,:));
end
fclose(fid);



