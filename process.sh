#!/bin/bash

#PBS -m be
#PBS -e ${t}.${n}n.${p}.${l}.${k}.error
#PBS -o ${t}.${n}n.${p}.${l}.${k}.output
#PBS -N matrix_power

num_nodes=$[$p*$l]
mpiexec -np $num_nodes ${HOME}/proj3/matrix_power -n $n -t $t -q $q -k $k -c -a
