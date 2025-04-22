#!/bin/bash -l
#PBS -N HYPRE_gpu_job
#PBS -l select=1
#PBS -l place=scatter
#PBS -l walltime=2:00:00
#PBS -l filesystems=home:eagle
#PBS -j oe
#PBS -q preemptable
#PBS -A CampSwift

# Modules:
module load PrgEnv-gnu
module load nvhpc-mixed
module load craype-accel-nvidia80
export MPICH_GPU_SUPPORT_ENABLED=1

# MPI and OpenMP settings
NNODES=1
NRANKS_PER_NODE=${NUM_OF_MESHES} # Meshes
NDEPTH=1
NTHREADS=1
export AFFINITY_RANKS_PER_GPU=${RANKS_PER_GPU}
export FDS_RANKS_PER_GPU=${RANKS_PER_GPU}

NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))
echo "NUM_OF_NODES= ${NNODES} TOTAL_NUM_RANKS= ${NTOTRANKS} RANKS_PER_NODE= ${NRANKS_PER_NODE} THREADS_PER_RANK= ${NTHREADS}"

FDSEXEC=~/firemodels/firex/fds/Build/ompi_gnu_linux/fds_ompi_gnu_linux
echo "Launching application..."
echo "       Input file: ${FDS_CASE_FILE}"
echo "    NUM_OF_MESHES: ${NUM_OF_MESHES}"
echo "          RS_CASE: ${RS_CASE}"
echo "        Exec file: ${FDSEXEC}"
echo "FDS_RANKS_PER_GPU: ${FDS_RANKS_PER_GPU}"
echo "        Directory: `pwd`"
echo "             Host: `hostname`"
echo "started running at `date`" >> ${FDS_BASENAME}.qlog
# For applications that internally handle binding MPI/OpenMP processes to GPUs
cd /home/isapchp/firemodels/firex/fds/Verification/GPU_Tests/HYPRE_GPU_SCALING
mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --depth=${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS=${NTHREADS} -env OMP_PLACES=threads ./set_affinity_gpu_polaris.sh ${FDSEXEC}  ${FDS_CASE_FILE}

