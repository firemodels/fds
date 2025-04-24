#!/bin/bash

#!/bin/bash
#SBATCH -A ENG146
#SBATCH -J test_gpu
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -t 1:00:00
#SBATCH -q debug
#SBATCH -N 1

module load PrgEnv-gnu
module load rocm
module load craype-accel-amd-gfx90a
export MPICH_GPU_SUPPORT_ENABLED=1
export FDS_RANKS_PER_GPU=8



FDSEXEC=~/firemodels/firex/fds/Build/ompi_gnu_linux/fds_ompi_gnu_linux
echo "Launching application..."
echo "       Input file: test_gpu.fds"
echo "        Exec file: ${FDSEXEC}"
echo "FDS_RANKS_PER_GPU: ${FDS_RANKS_PER_GPU}"
echo "        Directory: `pwd`"
echo "             Host: `hostname`"
echo "started running at `date`" >> test_gpu.qlog
cd /ccs/home/isapchp/firemodels/firex/fds/Verification/GPU_Tests
#srun -N 1 -n 8 -c 1 -m block:block --ntasks-per-gpu=8 --gpus-per-node=1 --gpu-bind=none  ${FDSEXEC} test_gpu.fds | sort
srun -n 8 -c 1 -m block:block  --gpus-per-task=1 --gpu-bind=closest  ${FDSEXEC} test_gpu.fds | sort
echo "finished running at `date`" >> test_gpu.qlog

