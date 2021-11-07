#!/bin/bash

export DB_PATH=./db
export OMP_NUM_THREADS=16

echo "./sample_ratio.exe ~/graph/m3.lg"
./sample_ratio.exe ~/graph/m3.lg

echo "./sample_ratio.exe ~/graph/com-orkut_r.lg"
./sample_ratio.exe ~/graph/com-orkut_r.lg

echo "./sample_ratio.exe ~/graph/uk-2005.lg"
./sample_ratio.exe ~/graph/uk-2005.lg
