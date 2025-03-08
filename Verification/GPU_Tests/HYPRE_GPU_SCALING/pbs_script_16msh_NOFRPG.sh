#!/bin/bash -l
#PBS -l select=1:system=polaris
#PBS -l place=scatter
#PBS -l walltime=0:30:00
#PBS -l filesystems=home:eagle
#PBS -j oe
#PBS -q preemptable
#PBS -A CampSwift

# Modules:
module load cudatoolkit-standalone/12.5.0 PrgEnv-gnu cray-libsci
# Enable GPU-MPI (if supported by application)
export MPICH_GPU_SUPPORT_ENABLED=1

# Change to working directory
cd $FIREMODELS/fds/Verification/GPU_Tests/HYPRE_GPU_SCALING

# MPI and OpenMP settings
NNODES=1
NRANKS_PER_NODE=16 # Meshes
NDEPTH=1
NTHREADS=1
export AFFINITY_RANKS_PER_GPU=4

NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))
echo "NUM_OF_NODES= ${NNODES} TOTAL_NUM_RANKS= ${NTOTRANKS} RANKS_PER_NODE= ${NRANKS_PER_NODE} THREADS_PER_RANK= ${NTHREADS}"

# For applications that need mpiexec to bind MPI ranks to GPUs
mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --depth=${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS=${NTHREADS} -env OMP_PLACES=threads ./set_affinity_gpu_polaris.sh $FIREMODELS/fds/Build/ompi_gnu_linux/fds_ompi_gnu_linux test_16mesh_NOFRPG.fds 
