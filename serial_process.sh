#!/bin/bash

#PBS -m be
#PBS -e ${t}.${n}n.${p}.${l}.error
#PBS -o ${t}.${n}n.${p}.${l}.output
#PBS -N matrix_power

num_nodes=$[$p*$l]
${HOME}/proj3/matrix_power_serial -n $n -t $t -q $q -c
