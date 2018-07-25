#!/bin/bash

echo "---------------------------------------------------------"
echo " Initializing maxwell relaxation test."
echo "---------------------------------------------------------"

RUN=./maxwell_relax
INPUT_FILE="./input_files/maxwell_relax.inp"

$RUN $INPUT_FILE
exit $?