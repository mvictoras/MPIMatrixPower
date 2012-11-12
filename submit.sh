#!/bin/bash

if ( ! getopts ":n:p:t:k:q:w:l:" opt); then
  echo "Usage: `basename $0` -n <number of numbers> -p <number of physical processors> -t <serial, mesh> -k <processors per node>"
  exit 1
fi
while getopts ":n:p:t:k:q:w:l:" opt; do
    case $opt in
    n)
      n=$OPTARG
      ;;
    p)
      p=$OPTARG
      ;;
    t)
      t=$OPTARG
      ;;
    k)
      k=$OPTARG
      ;;
    l)
      l=$OPTARG
      ;;
    w)
      w=$OPTARG
      ;;
    q)
      q=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

echo "Creating $n random numbers and using $p physical processors with $l processors per node $t $k $q"
qsub -V -l nodes=$p:ppn=$l -q student_long -v l=$l,n=$n,p=$p,t=$t,k=$k,q=$q process.sh 
