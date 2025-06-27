#!/bin/bash

BIN="./Binaires"

# Compilation

make all EXACTE=0 ECRITURE=0 ARRET=1

# Exécution

echo " "
echo "EXECUTIONS POUR PROBLEME-CHALEUR (mesure du temps)"
echo " "
sleep 1

$BIN/sequentiel-1 0 1 200 1 160000

for n in 1 2 4 6 8; do
    OMP_NUM_THREADS=$n $BIN/parallele-1 0 1 200 1 160000
done

for n in 1 2 4 6 8; do
    mpiexec -n $n $BIN/parallele-2 0 1 200 1 160000
done

for N in 200 400 600 800 1000 1200 1400; do
    $BIN/sequentiel-2 0 1 $N 1 $N
    $BIN/sequentiel-3 0 1 $N 1 $N
done