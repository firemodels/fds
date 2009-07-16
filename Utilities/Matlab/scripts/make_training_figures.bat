d:
cd \FDS-SMV\Utilities\Matlab
copy master_training_script.m temp.m
echo exit >> temp.m
matlab -r temp -nodesktop -nojvm -nosplash 
