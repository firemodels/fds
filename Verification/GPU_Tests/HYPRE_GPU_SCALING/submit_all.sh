#!/bin/bash

HPCSYS=$1

if [ "$HPCSYS" == "VISTA" ]; then
    export RANKS_PER_GPU=4
    export FDS_CASE_FILE=test_4mesh_RS4.fds
    sbatch -N 1 -n 4 VISTA_TACC_runjob.sh

    export RANKS_PER_GPU=2
    export FDS_CASE_FILE=test_8mesh_RS2.fds
    sbatch -N 4 -n 8 VISTA_TACC_runjob.sh

    export RANKS_PER_GPU=4
    export FDS_CASE_FILE=test_16mesh_RS4.fds
    sbatch -N 4 -n 16 VISTA_TACC_runjob.sh

    export RANKS_PER_GPU=8
    export FDS_CASE_FILE=test_32mesh_RS8.fds
    sbatch -N 4 -n 32 VISTA_TACC_runjob.sh

    export RANKS_PER_GPU=1
    export FDS_CASE_FILE=test_4mesh.fds
    sbatch -N 4 -n 4 VISTA_TACC_runjob.sh

    export RANKS_PER_GPU=1
    export FDS_CASE_FILE=test_8mesh_NOFRPG.fds
    sbatch -N 4 -n 8 VISTA_TACC_runjob.sh

    export RANKS_PER_GPU=1
    export FDS_CASE_FILE=test_16mesh_NOFRPG.fds
    sbatch -N 4 -n 16 VISTA_TACC_runjob.sh

    export RANKS_PER_GPU=1
    export FDS_CASE_FILE=test_32mesh_NOFRPG.fds
    sbatch -N 4 -n 32 VISTA_TACC_runjob.sh

elif [ "$HPCSYS" == "POLARIS" ]; then
    # All following job assume one node
    qsub -v RS_CASE=true,NUM_OF_MESHES=4,RANKS_PER_GPU=4,FDS_CASE_FILE=test_4mesh_RS4.fds POLARIS_ALCF_runjob.sh
    qsub -v RS_CASE=true,NUM_OF_MESHES=8,RANKS_PER_GPU=2,FDS_CASE_FILE=test_8mesh_RS2.fds POLARIS_ALCF_runjob.sh
    qsub -v RS_CASE=true,NUM_OF_MESHES=16,RANKS_PER_GPU=4,FDS_CASE_FILE=test_16mesh_RS4.fds POLARIS_ALCF_runjob.sh
    qsub -v RS_CASE=true,NUM_OF_MESHES=32,RANKS_PER_GPU=8,FDS_CASE_FILE=test_32mesh_RS8.fds POLARIS_ALCF_runjob.sh

    qsub -v NUM_OF_MESHES=4,RANKS_PER_GPU=1,FDS_CASE_FILE=test_4mesh.fds POLARIS_ALCF_runjob.sh
    qsub -v NUM_OF_MESHES=8,RANKS_PER_GPU=1,FDS_CASE_FILE=test_8mesh_NOFRPG.fds POLARIS_ALCF_runjob.sh
    qsub -v NUM_OF_MESHES=16,RANKS_PER_GPU=1,FDS_CASE_FILE=test_16mesh_NOFRPG.fds POLARIS_ALCF_runjob.sh
    qsub -v NUM_OF_MESHES=32,RANKS_PER_GPU=1,FDS_CASE_FILE=test_32mesh_NOFRPG.fds POLARIS_ALCF_runjob.sh

elif [ "$HPCSYS" == "FRONTIER" ]; then
    export NUM_OF_MESHES=4
    export RANKS_PER_GPU=4
    export FDS_CASE_FILE=test_4mesh_RS4.fds
    sbatch FRONTIER_ORNL_runjob.sh

    export NUM_OF_MESHES=8
    export RANKS_PER_GPU=2
    export FDS_CASE_FILE=test_8mesh_RS2.fds
    sbatch FRONTIER_ORNL_runjob.sh

    export NUM_OF_MESHES=16
    export RANKS_PER_GPU=4
    export FDS_CASE_FILE=test_16mesh_RS4.fds
    sbatch FRONTIER_ORNL_runjob.sh

    export NUM_OF_MESHES=32
    export RANKS_PER_GPU=8
    export FDS_CASE_FILE=test_32mesh_RS8.fds
    sbatch FRONTIER_ORNL_runjob.sh

    export NUM_OF_MESHES=4
    export RANKS_PER_GPU=1
    export FDS_CASE_FILE=test_4mesh.fds
    sbatch FRONTIER_ORNL_runjob.sh

    export NUM_OF_MESHES=8
    export RANKS_PER_GPU=1
    export FDS_CASE_FILE=test_8mesh_NOFRPG.fds
    sbatch FRONTIER_ORNL_runjob.sh

    export NUM_OF_MESHES=16
    export RANKS_PER_GPU=1
    export FDS_CASE_FILE=test_16mesh_NOFRPG.fds
    sbatch FRONTIER_ORNL_runjob.sh

    export NUM_OF_MESHES=32
    export RANKS_PER_GPU=1
    export FDS_CASE_FILE=test_32mesh_NOFRPG.fds
    sbatch FRONTIER_ORNL_runjob.sh

fi


