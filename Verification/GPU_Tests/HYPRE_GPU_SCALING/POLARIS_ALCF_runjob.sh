#!/bin/bash -l
#PBS -N test_gpu
#PBS -l select=1
#PBS -l place=scatter
#PBS -l walltime=0:30:00
#PBS -l filesystems=home:eagle
#PBS -j oe
#PBS -q debug-scaling
#PBS -A CampSwift

# Modules:
module load cudatoolkit-standalone/12.5.0 PrgEnv-gnu cray-libsci
# Enable GPU-MPI (if supported by application)
export MPICH_GPU_SUPPORT_ENABLED=1

# MPI and OpenMP settings
NNODES=1
NRANKS_PER_NODE=32 # Meshes
NDEPTH=1
NTHREADS=1
export AFFINITY_RANKS_PER_GPU=8
export FDS_RANKS_PER_GPU=${AFFINITY_RANKS_PER_GPU}

NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))
echo "NUM_OF_NODES= ${NNODES} TOTAL_NUM_RANKS= ${NTOTRANKS} RANKS_PER_NODE= ${NRANKS_PER_NODE} THREADS_PER_RANK= ${NTHREADS}"

FDSEXEC=~/firemodels/firex/fds/Build/ompi_gnu_linux/fds_ompi_gnu_linux
echo "Launching application..."
echo "       Input file: test_gpu.fds"
echo "        Exec file: ${FDSEXEC}"
echo "FDS_RANKS_PER_GPU: ${FDS_RANKS_PER_GPU}"
echo "        Directory: `pwd`"
echo "             Host: `hostname`"
echo "started running at `date`" >> test_gpu.qlog
# For applications that internally handle binding MPI/OpenMP processes to GPUs
cd /home/isapchp/firemodels/firex/fds/Verification/GPU_Tests/HYPRE_GPU_SCALING
mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --depth=${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS=${NTHREADS} -env OMP_PLACES=threads ./set_affinity_gpu_polaris.sh ${FDSEXEC} test_32mesh_RS8.fds

