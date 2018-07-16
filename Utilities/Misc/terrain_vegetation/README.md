## Instructions for use of define_vegetation.m script:

These instructions assume you have placed your *smv* and *fds* repos in the same directory. 

1. Go to *smv* repo and update it.

2. Navigate to the build directory of dem2fds: **smv/Build/dem2fds** and 
select your compilation target.
In burn, to compile with Intel ifort, select intel_linux_64 and execute
*./make_dem2fds.sh*
This will provide the executable *dem2fds_linux_64*.

3. Execute the shell script that will invoke *dem2fds*: 
In our case we are interested in the Gatlinburg 1000 m x 1000 m terrain case, 
defined in the *dem2fds* input script *Gatlinburg_1000m.in*. The script is
**Make_Gatlinburg_fds_input.sh**
This script runs dem2fds and produces both the OBST and GEOM version of
the terrain of interest and also image data in .png format.
For the Gatlinburg test case the name of the fds file using ZVALS (GEOM)
is Gatlinburg_1000m_g.fds

Note: Several examples on the execution of dem2fds can be found in:
**smv/Build/dem2fds/data**

4. Use the fds input file name (without the .fds extension) for geoms as the variable *casename*
on the **define_vegetation.m** Matlab script, and execute it.

The scripts tree terrain and size distribution parameters are:

- Total number of trees to distribute in domain: **N_TREES**

- Tree cone vertical distance from terrain provided by Zvals: **H_FROM_TERRAIN** (meters)

- Tree type : either 'block' or 'cone': **TREE_TYPE**

- Mean Tree cone base radius: **RADIUS** (meters)

- Maximum variation of cone base radius: **DRADIUS** (meters)

- Definition parameter for buffer layer from domain boundaries with no trees: **BUFF_DIST** (meters)

- Mean cone height: **HEIGHT** (meters)

- Maximum height variation: **DHEIGHT** (meters)

- Mass per volume for the trees: **MASS_PER_VOLUME**  (kg/m^3) 

- Particles per fds fluid cell: **N_PARTICLES_PER_CELL**

Note: If you rotate the final matlab plot you can see the tree
distribution in elevation.
