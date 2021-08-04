# Askervein Hill Input Files

This README documents the generation of FDS input files for the Askervein Hill series.

The [Askervein 1983](https://drive.google.com/file/d/1gb21m8G6irryDDNEoKsXwXBETwlGnccK/view?usp=sharing) gives the relevant information.

Thanks to Javier Sanz Rodrigo of the [Windbench](https://windbench.net/askervein-neutral) project for providing the topology data.  Details for generating the terrain data at 5 m resolution are given in the `firemodels/cad/Askervein_Hill/README.md` file.  The result of this process is the binary `Askervein_5m.gingeom` file that is used by the input files in this directory.

The grid is stretched (approx. 6%) in the vertical direction.  The transforms for this mapping are created using the Matlab script `fds/Utilities/Matlab/scripts/mesh_quad_stretch.m` with the following parameters:

Towers commented out for clarity:
The block of towers at the begining is commented out because they are not inside the domain of the mesh
HT Hill Top and CP Centre Point are commented out for clarity because they are not anemometer towers, they are locations
CP BSE 40 is commented out because it is the same tower as BSE 40.

Notes on tower lables as the relate to the table of towers in the report on page 144:
-towers labled t provide turbulence data, (DIR,UPWASH,SPEED,SIGMu,SIGMv,SIGMw,uv,uw,vw) and
-towers labled mf provide mean flow data (SPEED, DIRECTION, and sometimes FSUP(mean W velocity))
-ASW UK 30 m tower is a tower that provides vertical TU profile data, and it's tom height serves as ASW60
-SPEED CP FRG 17 m tower is a tower that provides vertical tu and MF profile data. IT is sometimes mislabled as 16m- and sometimes called CP' in the report.
-CP BSE 40 and BSE40 are the same tower.


The info below is likley no longer accurate, and will likley continue to change.

* 4 m resolution (vertical)
```
% Mesh size:
xs = -20.25; % Mesh block is defined from xs=1 to xf=xs+Lx=2.
Lx = 1000.0;

% Low side number m of cells with small uniform size DXS:
m  =     48;
DXS=     4.;

% Total number of cells:
Nx =     90;
```

* 8 m resolution
```
% Mesh size:
xs = -20.25; % Mesh block is defined from xs=1 to xf=xs+Lx=2.
Lx = 1000.0;

% Low side number m of cells with small uniform size DXS:
m  =     28;
DXS=     8.;

% Total number of cells:
Nx =     58;
```

* 16 m resolution (used for quick testing)
```
% Mesh size:
xs = -20.25; % Mesh block is defined from xs=1 to xf=xs+Lx=2.
Lx = 1000.0;

% Low side number m of cells with small uniform size DXS:
m  =     16;
DXS=     16.;

% Total number of cells:
Nx =     38;
```