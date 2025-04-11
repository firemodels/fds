#!/bin/bash
#SBATCH --job-name=ext_heartbeat_std_curve_python
#SBATCH --output=ext_heartbeat_std_curve_python.log
#SBATCH --error=ext_heartbeat_std_curve_python.err
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=99-99:99:99
export OMP_NUM_THREADS=1

# Load Python environment 
source $FIREMODELS/fds/.github/fds_python_env/bin/activate

echo "Running ext_heartbeat_std_curve_slurm.sh script..."
python ext_heartbeat_std_curve.py
echo "Done."
