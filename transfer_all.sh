#! /bin/sh

# transfer_all xyz transfiere todos los archivos necesarios al usuario
# xyz de atom

files='core.c core.h Makefile parameters.h tiny_md.c viz.c wtime.c wtime.h run.sh'

for f in $files; do
    scp $f cp2025_$1@atom.famaf.unc.edu.ar:/users/cp2025/cp2025_$1/tiny_md
done

