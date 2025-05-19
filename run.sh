#! /bin/sh

mod="" #! mod="perf stat"

if [ -z "$2" ] 
    then
        echo "Tenes que poner el valor de M y número de threads!!"
        echo "Por ej, con M=10:"
        echo "$> sh run.sh 10 2 4 8"
    else
        if [ ! -f log.csv ]
            then sh clear_log.sh > /dev/null
        fi
        M=$1
        shift
        make clean > /dev/null
        make M=$M tiny_md >/dev/null 2> /dev/null
        if [ ! -f tiny_md ]
            then    
                echo "Error de compilación (revisá los parámetros!)"
            else
                for nthreads in "$@"
                do
                    printf '%d,' "$nthreads" >> log.csv
                    OMP_NUM_THREADS=$nthreads $mod ./tiny_md
                    echo "------------------------------------------------------------------"
                done
        fi
fi