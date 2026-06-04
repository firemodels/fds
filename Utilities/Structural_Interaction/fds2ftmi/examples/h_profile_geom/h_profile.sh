mpiexec -n 1 ../../../../../Build/impi_intel_linux/fds_impi_intel_linux h_profile.fds 
%ANSYS% -b nolist -j h_profile_T_ansys -i read_geometry_h_profile.ans
python3 ../../source/pyfds2ftmi.py --chid h_profile --cell 0.05 --var 2 --avg mean:10 --code ansys --file h_profile_to_ansys.dat 
%ANSYS% -b nolist -j h_profile_T_ansys -i run_h_profile_T.ans
%ANSYS% -b nolist -j h_profile_M_ansys -i run_h_profile_M.ans 
%ANSYS% -b nolist -j h_profile_M_ansys -i plot_h_profile.ans 
export DISPLAY=localhost:0.0 #
display150 -d POSTSCRIPT -j h_profile_M_ansys #
epstopdf --outfile=h_profile_temp.pdf pscr000.eps  # 
epstopdf --outfile=h_profile_temp3d.pdf pscr001.eps  # 
epstopdf --outfile=h_profile_disp_15x.pdf pscr002.eps  # 
mv h_profile_temp.pdf ../../scripts/SCRIPT_FIGURES/h_profile_geom_temp.pdf #
mv h_profile_temp3d.pdf ../../scripts/SCRIPT_FIGURES/h_profile_geom_temp3d.pdf #
mv h_profile_disp_15x.pdf ../../scripts/SCRIPT_FIGURES/h_profile_geom_disp_15x.pdf #
mv "h_profile_to_ansys.dat" "h_profile_ftmi.txt" #
rm *.lock #
rm *.dat #
rm *_ansys* #
rm *.tmp #
rm *.eps #
