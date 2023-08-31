#!/bin/bash -l
#
# This pbs submission script runs test_gpu.fds defined with 8 meshes in Polaris.
# 1. Change the env var "FDSDIR" by the full path to your fds directory. 
# 2. Run the case typing:
#    $qsub pbs_script.sh
# 3. To see if it is running : $qstat -u YOUR_USER
#    Note the number of job (JOBID) in first column.
#    To check if gpus are engaged login to a node in use (you'll see them in pbs_script.sh.oJOBID) :
#    $ssh nodename
#    $nvidia-smi
#    To delete the job : $qdel JOBID
# 4. You can change the NTHREADS, and number of meshes in test_gpu.fds (for 4 meshes use select=1 node 
#    below, for 16 meshes select=4, etc.)
#
#PBS -l select=2:system=polaris
#PBS -l place=scatter
#PBS -l walltime=0:30:00
#PBS -l filesystems=home:eagle
#PBS -j oe
#PBS -q debug
#PBS -A CampSwift

# Modules:
module load cudatoolkit-standalone PrgEnv-gnu cray-libsci
# Enable GPU-MPI (if supported by application)
export MPICH_GPU_SUPPORT_ENABLED=1

# Where is FDS:
FDSDIR=/home/mnv/Firemodels_fork/fds

# Change to working directory
cd ${FDSDIR}/Verification/GPU_Tests

# MPI and OpenMP settings
NNODES=`wc -l < $PBS_NODEFILE`           # Defined from the select=2 nodes PBS line up top. Each node has 32 physical cores in Polaris.
NRANKS_PER_NODE=$(nvidia-smi -L | wc -l) # Picks up how many GPUs are present in the node. There are 4 in Polaris, so 4 MPI processes/node. 
NDEPTH=8                                 # Assign an MPI process every 8 physical cores (cores 0,7,15,23).
NTHREADS=8                               # 8 OpenMP threads per MPI process.

NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))
echo "NUM_OF_NODES= ${NNODES} TOTAL_NUM_RANKS= ${NTOTRANKS} RANKS_PER_NODE= ${NRANKS_PER_NODE} THREADS_PER_RANK= ${NTHREADS}"

# For applications that internally handle binding MPI/OpenMP processes to GPUs
mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --depth=${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS=${NTHREADS} -env OMP_PLACES=threads ${FDSDIR}/Build/mpich_gnu_polaris/fds_mpich_gnu_polaris test_gpu.fds -mat_type aijcusparse -vec_type cuda 

