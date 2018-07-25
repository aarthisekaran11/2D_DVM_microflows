#!/bin/bash

echo "---------------------------------------------------------"
echo " Initializing bkw relaxation test."
echo "---------------------------------------------------------"

RUN=./bkw_relax
INPUT_FILE="./input_files/bkw_relax.inp"

$RUN $INPUT_FILE
exit $?