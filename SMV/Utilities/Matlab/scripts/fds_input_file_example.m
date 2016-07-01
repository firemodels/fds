% McDermott
% 11 June 2015
% fds_input_file_example.m

close all
clear all

fid = fopen('chid.fds','wt');

tmp_str = ['&HEAD CHID=''chid'', TITLE=''whatever'' /']; fprintf(fid,'%s\n',tmp_str);

fprintf(fid,'%s\n','  '); % blank line

tmp_str = ['&MESH IJK=16,16,16, XB=0,1,0,1,0,1/']; fprintf(fid,'%s\n',tmp_str);

fprintf(fid,'%s\n','  '); % blank line

tmp_str = ['&TIME T_END=10./']; fprintf(fid,'%s\n',tmp_str);

fprintf(fid,'%s\n','  '); % blank line

tmp_str = ['&SURF ID=''fan'', VEL=-1., COLOR=''BLUE'' /']; fprintf(fid,'%s\n',tmp_str);

fprintf(fid,'%s\n','  '); % blank line

tmp_str = ['&VENT MB=''XMAX'', SURF_ID=''OPEN'' /']; fprintf(fid,'%s\n',tmp_str);
tmp_str = ['&VENT MB=''XMIN'', SURF_ID=''fan'' /'];  fprintf(fid,'%s\n',tmp_str);
tmp_str = ['&VENT MB=''YMAX'', SURF_ID=''OPEN'' /']; fprintf(fid,'%s\n',tmp_str);
tmp_str = ['&VENT MB=''YMIN'', SURF_ID=''OPEN'' /']; fprintf(fid,'%s\n',tmp_str);
tmp_str = ['&VENT MB=''ZMAX'', SURF_ID=''OPEN'' /']; fprintf(fid,'%s\n',tmp_str);

fprintf(fid,'%s\n','  '); % blank line

tmp_str = ['&SLCF PBY=0, QUANTITY=''VELOCITY'', VECTOR=.TRUE. /'];  fprintf(fid,'%s\n',tmp_str);

fprintf(fid,'%s\n','  '); % blank line

% write OBST lines

% create obst(s)
obst = struct('xb',[.1 .2 .1 .2 .1 .2]);

for i=1:length(obst)
    tmp_str = ['&OBST XB=',num2str(obst.xb(1)),',',num2str(obst.xb(2)),',',...
                           num2str(obst.xb(3)),',',num2str(obst.xb(4)),',',...
                           num2str(obst.xb(5)),',',num2str(obst.xb(6)),' /']; fprintf(fid,'%s\n',tmp_str);
end

fprintf(fid,'%s\n','  '); % blank line

tmp_str = ['&TAIL /']; fprintf(fid,'%s\n',tmp_str);

fclose(fid);



    