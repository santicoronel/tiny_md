#! /bin/sh

mod="" #! mod="perf stat"

if [ -z "$1" ]
    then
        echo "Poné el numero de threads!! (por ej $> sh run.sh 2 4 8)"
    else
        make clean && make tiny_md
        echo "N,tiempo[s],métrica(N^2/t)" > log.csv
        for nthreads in "$@"
        do
            OMP_NUM_THREADS=$nthreads $mod ./tiny_md
        done
fi