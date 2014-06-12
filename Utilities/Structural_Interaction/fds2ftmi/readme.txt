At this point the FTMI (Fire-Thermomechanical Interface) is running between FDS and ANSYS codes.
We DO NOT endorse the use of ANSYS or any other FEA code.

Procedures to evaluate FTMI:

1- Generate the FDS model;
	FTMI uses the boundary files, so FDS needs to write 2 boundary files for the interface geometry: 
	adiabatic surface temperature and heat transfer coefficient.

2- Genereate the mesh in ANSYS mesh generator;
	The geometry should be compatible with the OBSTs created for FDS simulation;
	It means same Cartesian coordinate system.

3- Export the mesh using the macro "ansys2ftmi.geom" or "ansys2ftmi_shell.geom"
	Maybe you need to open the macro file and edit it using the surface numbers where you want to prescribe the interface;
	The shell version creates interfaces at the both sides of shell structural members;
	This step will generate 2 files: nodes.dat and elements.dat;

4- Run fds2ftmi;
	To compile the code use the makefile provided with the source; 
	This step will generate 2 files per surface element and 1 output file;

5- Back to the ANSYS model, read te output from fds2ftmi ("output".dat);
	This step will generate the thermal input into ANSYS model.