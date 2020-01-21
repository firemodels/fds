function [WRITEFLG]=write_bingeom(GEOM,basedir,fdscase,geomname)

% UNIX fortran FDS.

IAXIS=1; JAXIS=2; KAXIS=3;
NOD1 =1; NOD2=2;  NOD3=3; 

intfmt = 'int32';   fourbyte = 4;
fltfmt = 'float64'; eightbyte= 8;

WRITEFLG = 1;

bingeom_dir_filename = [basedir fdscase '_' geomname '_OUT.bingeom'];

[fid]=fopen(bingeom_dir_filename,'wb');

% Load INT_ONE:
[tag]=fwrite(fid,1*fourbyte,intfmt);
[COUNT] = fwrite(fid,1,intfmt);
[tag]=fwrite(fid,1*fourbyte,intfmt);

% Load N_VERTS,N_FACES,N_SURF_ID,N_VOLUS:
[tag]=fwrite(fid,4*fourbyte,intfmt);
[COUNT]   = fwrite(fid,GEOM.N_VERTS,intfmt);
[COUNT]   = fwrite(fid,GEOM.N_FACES,intfmt);
[COUNT]   = fwrite(fid,GEOM.N_SURF_ID,intfmt);
[COUNT]   = fwrite(fid,GEOM.N_VOLUS,intfmt);
[tag]=fwrite(fid,4*fourbyte,intfmt);

% Load VERTS, FACES, SURFS, VOLUS:
[tag]=fwrite(fid,3*eightbyte*GEOM.N_VERTS,intfmt);
[COUNT]     = fwrite(fid,GEOM.VERTS,fltfmt);
if(COUNT ~= 3*GEOM.N_VERTS)
    disp(['Error writing vertices, 3*N_VERTS, Written = ' ...
          num2str(3*GEOM.N_VERTS) ', ' num2str(COUNT)])
end
[tag]=fwrite(fid,3*eightbyte*GEOM.N_VERTS,intfmt);
[tag]=fwrite(fid,3*fourbyte*GEOM.N_FACES,intfmt);
[COUNT]     = fwrite(fid,GEOM.FACES,intfmt);
if(COUNT ~= 3*GEOM.N_FACES)
    disp(['Error writing faces, 3*N_FACES, Written = ' ...
          num2str(3*GEOM.N_FACES) ', ' num2str(COUNT)])
end
[tag]=fwrite(fid,3*fourbyte*GEOM.N_FACES,intfmt);
[tag]=fwrite(fid,fourbyte*GEOM.N_FACES,intfmt);
[COUNT]     = fwrite(fid,GEOM.SURFS,intfmt);
if(COUNT ~= GEOM.N_FACES)
    disp(['Error writing SURF_ID per face, N_FACES, Written = ' ...
          num2str(GEOM.N_FACES) ', ' num2str(COUNT)])
end
[tag]=fwrite(fid,fourbyte*GEOM.N_FACES,intfmt);
if(GEOM.N_VOLUS > 0)
    [tag]=fwrite(fid,4*fourbyte*GEOM.N_VOLUS,intfmt);
    [COUNT]     = fwrite(fid,GEOM.VOLUS,intfmt);
    [tag]=fwrite(fid,4*fourbyte*GEOM.N_VOLUS,intfmt);
end
fclose(fid);

WRITEFLG=0;

return