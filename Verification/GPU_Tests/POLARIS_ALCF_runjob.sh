#!/bin/bash -l
#PBS -N test_gpu
#PBS -l select=1
#PBS -l place=scatter
#PBS -l walltime=0:30:00
#PBS -l filesystems=home:eagle
#PBS -j oe
#PBS -q debug
#PBS -A CampSwift

# Modules:
module load PrgEnv-gnu
module load nvhpc-mixed
module load craype-accel-nvidia80
export MPICH_GPU_SUPPORT_ENABLED=1

# MPI and OpenMP settings
NNODES=`wc -l < $PBS_NODEFILE`           # Defined from the select=2 nodes PBS line up top. Each node has 32 physical cores in Polaris.
NRANKS_PER_NODE=8 #$(nvidia-smi -L | wc -l) # Picks up how many GPUs are present in the node. There are 4 in Polaris, so 4 MPI processes/node. 
NDEPTH=4                                 # Assign an MPI process every 8 physical cores (cores 0,7,15,23).
NTHREADS=4                               # 8 OpenMP threads per MPI process.

NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))
echo "NUM_OF_NODES= ${NNODES} TOTAL_NUM_RANKS= ${NTOTRANKS} RANKS_PER_NODE= ${NRANKS_PER_NODE} THREADS_PER_RANK= ${NTHREADS}"

export FDS_RANKS_PER_GPU=8
FDSEXEC=~/firemodels/firex/fds/Build/ompi_gnu_linux/fds_ompi_gnu_linux
echo "Launching application..."
echo "       Input file: test_gpu.fds"
echo "        Exec file: ${FDSEXEC}"
echo "FDS_RANKS_PER_GPU: ${FDS_RANKS_PER_GPU}"
echo "        Directory: `pwd`"
echo "             Host: `hostname`"
echo "started running at `date`" >> test_gpu.qlog
# For applications that internally handle binding MPI/OpenMP processes to GPUs
cd /home/isapchp/firemodels/firex/fds/Verification/GPU_Tests
mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --depth=${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS=${NTHREADS} -env OMP_PLACES=threads ${FDSEXEC} test_gpu.fds 

