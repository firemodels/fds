#!/bin/bash -l
# This pbs submission script runs test_gpu.fds defined with 8 meshes.
# Run it with "sbatch VISTA_TACC_runjob.sh"
#SBATCH -J test_gpu_job          # Job name
#SBATCH -o test_gpu_job.%j.out   # Output file (%j = job ID)
#SBATCH -t 00:30:00              # Wall time
#SBATCH -N 1                     # Number of nodes
#SBATCH -n 8                     # Total tasks (4 per node)
#SBATCH -p gh-dev                # GPU partition (modify as needed)
#SBATCH -A CDA24014              # Project allocation

# Run your CUDA application
module load gcc/13.2.0
module load cuda
module load nvidia_math


# 1. Configure CUDA MPS environment
export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps
export CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-log

# Launch MPS daemon on all nodes
echo "Starting MPS daemon..."
export TACC_TASKS_PER_NODE=1     # Force 1 task/node
ibrun -np $SLURM_NNODES nvidia-cuda-mps-control -d
unset TACC_TASKS_PER_NODE
sleep 5                          # Wait for daemons to initialize

# 2. Run your CUDA application
export UCX_RCACHE_ENABLE=n
#export UCX_LOG_LEVEL=info
export OMP_NUM_THREADS=1
export FDS_RANKS_PER_GPU=8

echo "Launching application..."
echo "       Input file: test_gpu.fds"
echo "        Exec file: ${FDSEXEC}"
echo "FDS_RANKS_PER_GPU: ${FDS_RANKS_PER_GPU}"
echo "        Directory: `pwd`"
echo "             Host: `hostname`"
echo "started running at `date`" >> test_gpu.qlog
FDSEXEC=~/firemodels/firex/fds/Build/ompi_gnu_linux/fds_ompi_gnu_linux
ibrun ${FDSEXEC} test_gpu.fds   
echo "finished running at `date`" >> test_gpu.qlog

