#!/bin/bash
#SBATCH --job-name=ext_heartbeat_std_curve_python
#SBATCH --output=Controls/%x.log
#SBATCH --error=Controls/%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=99-99:99:99
export OMP_NUM_THREADS=1

# Load Python environment 
source $FIREMODELS/fds/.github/fds_python_env/bin/activate

# Set the working directory
cd $FIREMODELS/fds/Verification/Controls
echo "Running external heartbeat case..."
python ext_heartbeat_std_curve.py
echo "Done."
