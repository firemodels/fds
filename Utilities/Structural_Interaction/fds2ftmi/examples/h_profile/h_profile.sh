../../../../../Build/impi_intel_linux_64/fds_impi_intel_linux_64 h_profile.fds 
%ANSYS% -b nolist -j h_profile_T_ansys -i read_geometry_h_profile.ans
../../intel_linux_64/fds2ftmi_linux_64 h_profile 0.05 2 0 30 2 10 0 1 h_profile_to_ansys 
%ANSYS% -b nolist -j h_profile_T_ansys -i run_h_profile_T.ans
%ANSYS% -b nolist -j h_profile_M_ansys -i run_h_profile_M.ans 
%ANSYS% -b nolist -j h_profile_M_ansys -i plot_h_profile.ans 
export DISPLAY=localhost:0.0 #
display150 -d POSTSCRIPT -j h_profile_M_ansys #
epstopdf --outfile=h_profile_temp.pdf pscr000.eps  # 
epstopdf --outfile=h_profile_temp3d.pdf pscr001.eps  # 
epstopdf --outfile=h_profile_disp_15x.pdf pscr002.eps  # 
mv h_profile_temp.pdf ../../scripts/SCRIPT_FIGURES/h_profile_temp.pdf #
mv h_profile_temp3d.pdf ../../scripts/SCRIPT_FIGURES/h_profile_temp3d.pdf #
mv h_profile_disp_15x.pdf ../../scripts/SCRIPT_FIGURES/h_profile_disp_15x.pdf #
rm *.lock #
rm *.dat #
rm *_ansys* #
rm *.tmp #
rm *.eps #
