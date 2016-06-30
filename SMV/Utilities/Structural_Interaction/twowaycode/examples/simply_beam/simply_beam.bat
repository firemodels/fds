%ANSYS% -b nolist -j simply_beam_ansys -i read_geometry_simply_beam.ans -o simply_beam
python simply_beam.py
%ANSYS% -b -j simply_beam_step1 -i simply_beam_step1.ans -o simply_beam
%ANSYS% -b -j simply_beam_step2 -i simply_beam_step2.ans -o simply_beam
..\..\intel_win_64\twowaycode_win_64 input_simply_beam_win.2way 
%ANSYS% -b -j simply_beam_ansys_mechanical -i export_data_simply_beam.ans -o simply_beam
%ANSYS% -b -j simply_beam_ansys_mechanical -i plot_simply_beam.ans -o simply_beam 
set DISPLAY=localhost:0.0 
displayw -d POSTSCRIPT -j simply_beam_ansys 
epstopdf --outfile=simply_beam_temp.pdf pscr000.eps  
move simply_beam_temp.pdf ..\..\scripts\SCRIPT_FIGURES\simply_beam_temp.pdf 
del *.lock 
del *.dat 
del *_ansys* 
del *.tmp 
del *.eps 
del *.log 
del *.err 
del *.db 
del *.dbb 