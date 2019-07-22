#!/bin/bash

SCRIPT="$1"
shift
source ~/.bashrc
"$SCRIPT" $@
