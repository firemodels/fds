function [READFLG,GEOM]=load_bingeom(basedir,fdscase,geomname)

% UNIX fortran FDS.

IAXIS=1; JAXIS=2; KAXIS=3;
NOD1 =1; NOD2=2;  NOD3=3; 

intfmt = 'int32';
fltfmt = 'float64';

READFLG = 1;

bingeom_dir_filename = [basedir fdscase '_' geomname '.bingeom'];

[fid]=fopen(bingeom_dir_filename,'rb');

% Load INT_ONE:
[tag]=fread(fid,1,intfmt);

[INT_ONE] = fread(fid,1,intfmt);
 
[tag]=fread(fid,1,intfmt);

% Load N_VERTS,N_FACES,N_SURF_ID,N_VOLUS:
[tag]=fread(fid,1,intfmt);
[GEOM.N_VERTS]   = fread(fid,1,intfmt);
[GEOM.N_FACES]   = fread(fid,1,intfmt);
[GEOM.N_SURF_ID] = fread(fid,1,intfmt);
[GEOM.N_VOLUS]   = fread(fid,1,intfmt);
[tag]=fread(fid,1,intfmt);

% Load VERTS, FACES, SURFS, VOLUS:
[tag]=fread(fid,1,intfmt);
[GEOM.VERTS]     = fread(fid,3*GEOM.N_VERTS,fltfmt);
[tag]=fread(fid,1,intfmt);
[tag]=fread(fid,1,intfmt);
[GEOM.FACES]     = fread(fid,3*GEOM.N_FACES,intfmt);
[tag]=fread(fid,1,intfmt);
[tag]=fread(fid,1,intfmt);
[GEOM.SURFS]     = fread(fid,1*GEOM.N_FACES,intfmt);
[tag]=fread(fid,1,intfmt);
if(GEOM.N_VOLUS > 0)
    [tag]=fread(fid,1,intfmt);
    [GEOM.VOLUS]     = fread(fid,4*GEOM.N_VOLUS,intfmt);
    [tag]=fread(fid,1,intfmt);
end
fclose(fid);

% Build XYZ and WSELEM:
GEOM.XYZ=zeros(GEOM.N_VERTS,KAXIS);
for IVERT=1:GEOM.N_VERTS
    GEOM.XYZ(IVERT,IAXIS:KAXIS) = GEOM.VERTS(KAXIS*(IVERT-1)+1:KAXIS*IVERT);
end
GEOM.WSELEM=zeros(GEOM.N_FACES,NOD3);
for IFACE=1:GEOM.N_FACES
    GEOM.WSELEM(IFACE,NOD1:NOD3)= GEOM.FACES(NOD3*(IFACE-1)+1:NOD3*IFACE);
end

READFLG=0;

return