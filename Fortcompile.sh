#! /bin/bash

ifort -O3 ./src/preprocess/*.f90 -o ./preproces
echo '***compile successfully -> preproces*^_^'
ifort -O3 -heap-arrays ./src/solver/*.f90 -o ./solver
echo '***compile successfully -> solver*^_^'
ifort -O3 ./src/postprocess/*.f90 -o ../postproces
echo '***compile successfully -> postproces*^_^'