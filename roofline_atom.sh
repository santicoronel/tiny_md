#! /bin/sh

tiny_md_path=tiny_md ; #aca va tu path, yo tengo todo en tiny_md/

pwd=$PWD ;
cd $tiny_md_path ;
make clean ;
make M=$1 tiny_md;
cd $pwd ;

mkdir roofline ;

source /opt/intel/oneapi/setvars.sh ;
export OMP_NUM_THREADS=$2;
advixe-cl -collect roofline -project-dir ./roofline -- $tiny_md_path/tiny_md ;
advixe-cl -report roofline -project-dir ./roofline --report-output ./roofline.html ;

rm -rf roofline ;