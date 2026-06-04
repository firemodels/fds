mpiexec -n 1 ../../../../../Build/impi_intel_linux/fds_impi_intel_linux simple_panel_fire.fds 
%ANSYS% -b nolist -j simple_panel_fire_ansys -i read_geometry_fire.ans
python3 ../../source/pyfds2ftmi.py --chid simple_panel_fire --cell 0.05 --var 2 --code ansys --file simple_panel_fire_to_ansys.dat 
%ANSYS% -b nolist -j simple_panel_fire_ansys -i run_simple_panel_fire.ans 
%ANSYS% -b nolist -j simple_panel_fire_ansys -i plot_fire.ans 
export DISPLAY=localhost:0.0 #
display150 -d POSTSCRIPT -j simple_panel_fire_ansys #
epstopdf --outfile=simple_panel_fire_ansys.pdf pscr000.eps  # 
mv simple_panel_fire_ansys.pdf ../../scripts/SCRIPT_FIGURES/simple_panel_fire_ansys.pdf #
rm *.lock #
rm *.dat #
rm *_ansys* #
rm *.tmp #
rm *.eps #
