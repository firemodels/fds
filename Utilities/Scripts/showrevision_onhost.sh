#!/bin/bash

directory=$1
host=$2

ssh -q $host \( cd \~/$directory \; git describe --dirty --long  \)
