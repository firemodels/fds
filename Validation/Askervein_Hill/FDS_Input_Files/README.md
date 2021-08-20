# Askervein Hill Input Files

This README documents the generation of FDS input files for the Askervein Hill series.

The [Askervein 1983](https://drive.google.com/file/d/1gb21m8G6irryDDNEoKsXwXBETwlGnccK/view?usp=sharing) gives the relevant information.

Thanks to Javier Sanz Rodrigo of the [Windbench](https://windbench.net/askervein-neutral) project for providing the topology data.  Details for generating the terrain data at 5 m resolution are given in the `firemodels/cad/Askervein_Hill/README.md` file.  The result of this process is the binary `Askervein_5m.gingeom` file that is used by the input files in this directory.

The grid is stretched (approx. 6%) in the vertical direction.  The transforms for this mapping are created using the Matlab script `fds/Utilities/Matlab/scripts/mesh_quad_stretch.m` with the following parameters:

Towers commented out for clarity:
The block of towers at the begining is commented out because they are not inside the domain of the mesh
HT Hill Top and CP Centre Point are commented out for clarity because they are not anemometer towers, they are locations

Notes on tower lables as the relate to the table of towers in the report on page 144:
-towers labled t provide turbulence data, (DIR,UPWASH,SPEED,SIGMu,SIGMv,SIGMw,uv,uw,vw) and
-towers labled mf provide mean flow data (SPEED, DIRECTION, and sometimes FSUP(mean W velocity))
-ASW UK 30 m tower is a tower that provides vertical TU profile data, and it's tom height serves as ASW60
-SPEED CP FRG 17 m tower is a tower that provides vertical tu and MF profile data. 
It is actually 16m (it is labled 17m in the table of towers and in all of FDS) and is often called CP' in the report.
(The acutal FRG 17m tower is at RS and is a mean flow tower with it's highest cup anemometer and tempreture sensor at 16.9m)
-CP BSE 40 and BSE40 are the same tower.

Notes on the confidence and reasoning of Identification of remaining towers and data as they are used in the Askervein Hill Matlab preprocesing script:
-At the HT location, there is only one 10m tower able output tu data, so it should be the device recording data for the 0 position in the ASW-ANE line 
-At the HT location, there is only one 10m tower noted output mf data, so it should be the device recording data for the 0 position in the BSE-BNW line 

-At the CP location, there are many differnt graphs and tables and a few different towers that refrence data from this location. 
The AASW-AANE and BNW-BSE mean flow graphs have two datapoints at this location(likely the same two devices,) The CP' data  in Table A1.4 on page 87 
should be from the FRG17m tower, and the data labled CP MF (not CP') on this page should be from the tower labled CP FRG mf.
Because of this, I beleive that the two datapoints at CP in the Mean flow Graphs are the BSE40 and UK mf towers, and the script proceses it as such.

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