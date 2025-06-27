#!/bin/bash

BIN="./Binaires"

# Compilation

make all EXACTE=0 ECRITURE=0 ARRET=1

# Exécution

echo " "
echo "EXECUTIONS POUR PROBLEME-ONDES (mesure du temps)"
echo " "
sleep 1


for N in 5000 10000 20000; do

    $BIN/sequentiel-1 0 1 $N 1 $N

done