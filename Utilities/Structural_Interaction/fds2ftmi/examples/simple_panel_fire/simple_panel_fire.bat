..\..\..\..\..\Build\impi_intel_win_openmp\fds_impi_intel_win_openmp simple_panel_fire.fds 
%ANSYS% -b nolist -j simple_panel_fire_ansys -i read_geometry_fire.ans -o simple_panel_fire_ansys 
python ..\..\source\pyfds2ftmi.py --chid simple_panel_fire --cell 0.05 --var 2 --code ansys --file simple_panel_fire_to_ansys.dat 
%ANSYS% -b nolist -j simple_panel_fire_ansys -i run_simple_panel_fire.ans -o simple_panel_fire_ansys  
%ANSYS% -b nolist -j simple_panel_fire_ansys -i plot_fire.ans -o simple_panel_fire_ansys  
set DISPLAY=localhost:0.0 
displayw -d POSTSCRIPT -j simple_panel_fire_ansys 
epstopdf --outfile=simple_panel_fire_ansys.pdf pscr000.eps  
move simple_panel_fire_ansys.pdf ..\..\scripts\SCRIPT_FIGURES\simple_panel_fire_ansys.pdf 
del *.loc 
del *.dat 
del *_ansys* 
del *.tmp 
del *.eps 