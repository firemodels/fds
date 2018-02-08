#!/bin/bash

directory=$1
host=$2

cd ~/$directory
git describe --dirty --long
