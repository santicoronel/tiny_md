#! /bin/sh

# transfer xyz file1 file2 ... filen copia todos los archivos filei en el usuario
# xyz de atom

usr=$1 ;
shift ;
scp $@ cp2025_$usr@atom.famaf.unc.edu.ar:/users/cp2025/cp2025_$usr