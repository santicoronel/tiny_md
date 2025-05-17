#! /bin/sh
if [ -z "$1" ]
    then
        echo "Pone el numero de threads!! (por ej $> sh run.sh 2)"
    else 
        make clean && make tiny_md && OMP_NUM_THREADS=$1 ./tiny_md
fi