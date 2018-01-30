#!/bin/bash
export NOPAUSE=1
args=$0
DIR=$(dirname "${args}")
$DIR/make_bundle.sh
