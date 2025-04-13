#!/bin/bash

HPCSYS=$1

if [ "$HPCSYS" == "VISTA" ]; then
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
    sbatch -N 1 -n 4 VISTA_TACC_runjob.sh

    export RANKS_PER_GPU=1
    export FDS_CASE_FILE=test_8mesh_NOFRPG.fds
    sbatch -N 1 -n 8 VISTA_TACC_runjob.sh

    export RANKS_PER_GPU=1
    export FDS_CASE_FILE=test_16mesh_NOFRPG.fds
    sbatch -N 1 -n 16 VISTA_TACC_runjob.sh

    export RANKS_PER_GPU=1
    export FDS_CASE_FILE=test_32mesh_NOFRPG.fds
    sbatch -N 1 -n 32 VISTA_TACC_runjob.sh

elif [ "$HPCSYS" == "POLARIS" ]; then
    qsub pbs_script_4msh.sh
    qsub pbs_script_8msh_RS2.sh
    qsub pbs_script_8msh_NOFRPG.sh
    qsub pbs_script_16msh_RS4.sh
    qsub pbs_script_16msh_NOFRPG.sh
    qsub pbs_script_32msh_RS8.sh
    qsub pbs_script_32msh_NOiFRPG.sh
fi


