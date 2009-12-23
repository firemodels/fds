cd Controls
smokeview -runscript activate_vents
cd ..
cd Detectors
smokeview -runscript beam_detector
cd ..
cd Fires
smokeview -runscript room_fire
cd ..
cd Flowfields
smokeview -runscript helium_2d
smokeview -runscript sawtooth
cd ..
cd Miscellaneous
smokeview -runscript pyramid
cd ..
cd NS_Analytical_Solution
smokeview -runscript ns2d_64
cd ..
cd Pressure_Effects
smokeview -runscript pressure_boundary
cd ..
cd Sprinklers_and_Sprays
smokeview -runscript cascade

Rem generate Smokeview figures
cd ..
call Make_SMV_Pictures


