#!/bin/bash

#!/bin/bash
#SBATCH -A ENG146
#SBATCH -J test_gpu
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -t 16:00:00
#SBATCH -p extended
#SBATCH -N 1

module load PrgEnv-gnu
module load rocm
module load craype-accel-amd-gfx90a
export MPICH_GPU_SUPPORT_ENABLED=1
export FDS_RANKS_PER_GPU=8


# MPI and OpenMP settings
NNODES=1
NRANKS_PER_NODE=${NUM_OF_MESHES} # Meshes
NTHREADS=1

NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))
echo "NUM_OF_NODES= ${NNODES} TOTAL_NUM_RANKS= ${NTOTRANKS} RANKS_PER_NODE= ${NRANKS_PER_NODE} THREADS_PER_RANK= ${NTHREADS}"

export FDS_RANKS_PER_GPU=${RANKS_PER_GPU}


if [ "$RANKS_PER_GPU" -eq 1 ]; then
  CPARAM=1
  NTASKS_PER_GPU=7
else
  CPARAM=$(( 8 / RANKS_PER_GPU ))
  if [ "$CPARAM" -lt 1 ]; then
    CPARAM=1
  elif [ "$CPARAM" -gt 7 ]; then
    CPARAM=7
  fi
  NTASKS_PER_GPU=${RANKS_PER_GPU}
  if [ "$NTASKS_PER_GPU" -gt 7 ]; then
    NTASKS_PER_GPU=7
  fi
fi

FDSEXEC=~/firemodels/firex/fds/Build/ompi_gnu_linux/fds_ompi_gnu_linux
FDS_BASENAME="${FDS_CASE_FILE%.fds}"  # removes .fds suffix
echo "Launching application..."
echo "       Input file: ${FDS_CASE_FILE}"
echo "    NUM_OF_MESHES: ${NUM_OF_MESHES}"
echo "        Exec file: ${FDSEXEC}"
echo "FDS_RANKS_PER_GPU: ${FDS_RANKS_PER_GPU}"
echo "           CPARAM: ${CPARAM}"
echo "   NTASKS_PER_GPU: ${NTASKS_PER_GPU}"
echo "        Directory: `pwd`"
echo "             Host: `hostname`"
echo "started running at `date`" >> ${FDS_BASENAME}.qlog
cd /ccs/home/isapchp/firemodels/firex/fds/Verification/GPU_Tests/HYPRE_GPU_SCALING
if [ "$NTASKS_PER_GPU" -eq 7 ]; then
srun -n ${NTOTRANKS} -c ${CPARAM} -m block:block  --gpu-bind=closest  ${FDSEXEC} ${FDS_CASE_FILE}
else
srun -n ${NTOTRANKS} -c ${CPARAM} -m block:block --ntasks-per-gpu=${NTASKS_PER_GPU} --gpu-bind=closest  ${FDSEXEC} ${FDS_CASE_FILE}
fi
echo "finished running at `date`" >> ${FDS_BASENAME}.qlog

