#!/bin/bash

isRS=${RS_CASE}
num_gpus=4

# BEFORE USING THIS SCRIPT - DEFINE FDS_RANKS_PER_GPU in your submission script:
# export FDS_RANKS_PER_GPU=X
# where X is such that NRANKS_PER_NODE / X is not greater than 4.

# Need to assign GPUs in reverse order due to topology.
# See Polaris Device Affinity Information:
# https://www.alcf.anl.gov/support/user-guides/polaris/hardware-overview/machine-overview/index.html

# Packed in Rank:
if [[ "$isRS" == "true" ]]; then
    gpu=$((num_gpus - 1 - PMI_LOCAL_RANK / AFFINITY_RANKS_PER_GPU))
else
    gpu=$((num_gpus - 1 - PMI_LOCAL_RANK % num_gpus))
fi
export CUDA_VISIBLE_DEVICES=$gpu
echo "RANK=${PMI_RANK} LOCAL_RANK=${PMI_LOCAL_RANK} gpu=${gpu}"
exec "$@"


