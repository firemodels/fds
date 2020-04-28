% McDermott
% 12-6-2019
% parse_fds_out.m
%
% The basic plan is to start building this function and adding
% new values as needed.  Follow the examples below.
%
% Note that if someone has changed even one space in the formating of the
% .out file, then there will be problems with this parser.  So, beware
% that it is fragile.  Look at your .out file and fix the strcmp lines
% as needed.
%
% Usage: Vals = parse_fds_out([outdir,'chid.out'],'time step wall clock time','heat of combustion')

function Vals = parse_fds_out(filename,varargin)

Vals = zeros(1,length(varargin));

I_Time_Step_Wall_Clock_Time     = find(strcmp(varargin,'time step wall clock time'));
I_Total_Elapsed_Wall_Clock_Time = find(strcmp(varargin,'total elapsed wall clock time'));
I_Stoich_Coef                   = find(strcmp(varargin,'stoich coef'));
I_Heat_of_Combustion            = find(strcmp(varargin,'heat of combustion'));

fid = fopen(filename,'r+');

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end

    if I_Time_Step_Wall_Clock_Time>0
        if strfind(tline,'Time Stepping Wall Clock Time')
            % disp(tline)
            C = strsplit(tline);
            Vals(I_Time_Step_Wall_Clock_Time)=abs(str2num(C{end}));
        end
    end

    if I_Total_Elapsed_Wall_Clock_Time>0
        if strfind(tline,'Total Elapsed Wall Clock Time')
            % disp(tline)
            C = strsplit(tline);
            Vals(I_Total_Elapsed_Wall_Clock_Time)=abs(str2num(C{end}));
        end
    end

    if I_Stoich_Coef>0
        if strfind(tline,'Tracked (Lumped) Species Stoich. Coeff.')
            %disp(tline)
            tline = fgetl(fid); % dummy read
            tline = fgetl(fid); % read the AIR line
            C = strsplit(tline);
            Vals(I_Stoich_Coef)=abs(str2num(C{end}));
        end
    end

    if I_Heat_of_Combustion>0
        if strfind(tline,'Heat of Combustion (kJ/kg)')
            % disp(tline)
            tline = fgetl(fid); % read the FUEL line
            C = strsplit(tline);
            Vals(I_Heat_of_Combustion)=abs(str2num(C{end}));
        end
    end

end
fclose(fid);

return

