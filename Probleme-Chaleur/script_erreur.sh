#!/bin/bash

BIN="./Binaires"

# Compilation

make all EXACTE=1 ECRITURE=0 ARRET=0

# Exécution

echo " "
echo "EXECUTIONS POUR PROBLEME-CHALEUR (calcul des erreurs)"
echo " "
sleep 1

for h in 0.1 0.05 0.02 0.01; do
    for h_t in 0.0001 0.00005 0.00002 0.00001; do

        $BIN/sequentiel-1 1 1 $h 1 $h_t

        for n in 1 2 4 6 8; do
            OMP_NUM_THREADS=$n $BIN/parallele-1 1 1 $h 1 $h_t
        done

        for n in 1 2 4 6 8; do
            mpiexec -n $n $BIN/parallele-2 1 1 $h 1 $h_t
        done

        $BIN/sequentiel-2 1 1 $h 1 $h_t

        $BIN/sequentiel-3 1 1 $h 1 $h_t

    done
done

for h in 0.1 0.05 0.02 0.01; do
    for h_t in 0.1 0.05 0.02 0.01; do

        $BIN/sequentiel-2 1 1 $h 1 $h_t

        $BIN/sequentiel-3 1 1 $h 1 $h_t

    done
done