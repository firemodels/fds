mpiexec -n 1 ..\..\..\..\..\Build\impi_intel_win_openmp\fds_impi_intel_win_openmp h_profile.fds 
%ANSYS% -b nolist -j h_profile_T_ansys -i read_geometry_h_profile.ans -o h_profile_ansys 
python ..\..\source\pyfds2ftmi.py --chid h_profile --cell 0.05 --var 2 --avg mean:10 --code ansys --file h_profile_to_ansys.dat 
%ANSYS% -b nolist -j h_profile_T_ansys -i run_h_profile_T.ans -o h_profile_ansys 
%ANSYS% -b nolist -j h_profile_M_ansys -i run_h_profile_M.ans -o h_profile_ansys_1 
%ANSYS% -b nolist -j h_profile_M_ansys -i plot_h_profile.ans -o h_profile_ansys_2
set DISPLAY=localhost:0.0 
displayw -d POSTSCRIPT -j h_profile_M_ansys 
epstopdf --outfile=h_profile_temp.pdf pscr000.eps  
epstopdf --outfile=h_profile_temp3d.pdf pscr001.eps  
epstopdf --outfile=h_profile_disp_15x.pdf pscr002.eps  
move h_profile_temp.pdf ..\..\scripts\SCRIPT_FIGURES\h_profile_geom_temp.pdf 
move h_profile_temp3d.pdf ..\..\scripts\SCRIPT_FIGURES\h_profile_geom_temp3d.pdf 
move h_profile_disp_15x.pdf ..\..\scripts\SCRIPT_FIGURES\h_profile_geom_disp_15x.pdf 
ren "h_profile_to_ansys.dat" "h_profile_ftmi.txt"
del *.lock 
del *.dat 
del *_ansys* 
del *.tmp 
del *.eps 