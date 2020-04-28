..\..\..\..\..\Build\impi_intel_win_64\fds_impi_win_64 h_profile.fds 
%ANSYS% -b nolist -j h_profile_T_ansys -i read_geometry_h_profile.ans -o h_profile_ansys 
..\..\intel_win_64\fds2ftmi_win_64 h_profile 0.05 2 0 30 2 10 0 1 h_profile_to_ansys 
%ANSYS% -b nolist -j h_profile_T_ansys -i run_h_profile_T.ans -o h_profile_ansys 
%ANSYS% -b nolist -j h_profile_M_ansys -i run_h_profile_M.ans -o h_profile_ansys 
%ANSYS% -b nolist -j h_profile_M_ansys -i plot_h_profile.ans -o h_profile_ansys
set DISPLAY=localhost:0.0 
displayw -d POSTSCRIPT -j h_profile_M_ansys 
epstopdf --outfile=h_profile_temp.pdf pscr000.eps  
epstopdf --outfile=h_profile_temp3d.pdf pscr001.eps  
epstopdf --outfile=h_profile_disp_15x.pdf pscr002.eps  
move h_profile_temp.pdf ..\..\scripts\SCRIPT_FIGURES\h_profile_temp.pdf 
move h_profile_temp3d.pdf ..\..\scripts\SCRIPT_FIGURES\h_profile_temp3d.pdf 
move h_profile_disp_15x.pdf ..\..\scripts\SCRIPT_FIGURES\h_profile_disp_15x.pdf 
del *.lock 
del *.dat 
del *_ansys* 
del *.tmp 
del *.eps 