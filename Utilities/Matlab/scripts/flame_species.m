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

primitive_struct=importdata([FDS_Output_Files,'methane_flame_primitive_devc.csv']);
primitive_data=primitive_struct.data;
primitive_text=primitive_struct.textdata;

lumped_struct=importdata([FDS_Output_Files,'methane_flame_lumped_devc.csv']);
lumped_data=lumped_struct.data;
lumped_text=lumped_struct.textdata;

header1=cat(2,primitive_text(1),lumped_text(1,2:4));
header2=cat(2,primitive_text(2,:),lumped_text(2,2:4));
data=cat(2,primitive_data,lumped_data(:,2:4));

%-----------------------
% Write new files
%-----------------------

fid = fopen([FDS_Output_Files,'methane_flame_lumpedprimitive.csv'],'wt','n');

fprintf(fid,'%s, %s, %s, %s, %s, %s, %s, %s\n',header1{1,:});
fprintf(fid,'\n%s, %s, %s, %s, %s, %s, %s, %s\n',header2{1,:});
for i=1:length(primitive_data)
    fprintf(fid,'\n%f, %f, %f, %f, %f, %f, %f, %f\n',data(i,:));
end
fclose(fid);