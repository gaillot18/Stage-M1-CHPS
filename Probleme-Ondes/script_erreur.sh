#!/bin/bash

BIN="./Binaires"

# Compilation

make all EXACTE=1 ECRITURE=0 ARRET=0

# Exécution

echo " "
echo "EXECUTIONS POUR PROBLEME-ONDES (calcul des erreurs)"
echo " "
sleep 1

for h in 0.01 0.005 0.002 0.001; do
    for h_t in 0.01 0.005 0.002 0.001; do

        $BIN/sequentiel-1 1 1 $h 1 $h_t

    done
done