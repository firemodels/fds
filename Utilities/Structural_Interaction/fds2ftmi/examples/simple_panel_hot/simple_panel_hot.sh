../../../../../Build/impi_intel_linux_64/fds_impi_intel_linux_64 simple_panel_hot.fds 
%ANSYS% -b nolist -j simple_panel_hot_ansys -i read_geometry_hot.ans
../../intel_linux_64/fds2ftmi_linux_64 simple_panel_hot 0.05 2 0 600 2 0 0 1 simple_panel_hot_to_ansys 
%ANSYS% -b nolist -j simple_panel_hot_ansys -i run_simple_panel_hot.ans 
%ANSYS% -b nolist -j simple_panel_hot_ansys -i plot_hot.ans 
export DISPLAY=localhost:0.0 #
display150 -d POSTSCRIPT -j simple_panel_hot_ansys #
epstopdf --outfile=simple_panel_hot_ansys.pdf pscr000.eps  # 
mv simple_panel_hot_ansys.pdf ../../scripts/SCRIPT_FIGURES/simple_panel_hot_ansys.pdf #
rm *.lock #
rm *.dat #
rm *_ansys* #
rm *.tmp #
rm *.eps #
