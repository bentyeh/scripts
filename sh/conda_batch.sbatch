#!/bin/bash

ENV="$1"
SCRIPT="$2"
shift 2

source ~/.bashrc
conda activate $ENV
"$SCRIPT" $@
