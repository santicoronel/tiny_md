#! /bin/sh

# roofline 007 10 8 corre en el usuario 007 de atom, con M=10 y 8 cores
# configura tu path en roofline_atom.sh!!

ssh cp2025_$1@atom.famaf.unc.edu.ar 'bash -s' < roofline_atom.sh $2 $3
scp cp2025_$1@atom.famaf.unc.edu.ar:roofline.html .