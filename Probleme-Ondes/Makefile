CC = clang
CCP = mpicc

BIN = ./Binaires
LIB = ./Librairies
OBJ = ./Objets
SRC = ./Sources
TXT = ./Textes

FCC = Fonctions-communes
BS = base
PRL1 = parallele-1
PRL2 = parallele-2
SQT1 = sequentiel-1
SQT2 = sequentiel-2

FLAGS = -g0 -Wall -O3
MPI_FLAGS = -DUSE_MPI
OMP_FLAGS = -fopenmp
ifneq ($(EXACTE),0)
  FLAGS += -DEXACTE
endif
ifneq ($(ECRITURE),0)
  FLAGS += -DECRITURE
endif
ifneq ($(ARRET),0)
  FLAGS += -DARRET
endif

.SUFFIXES: .o .c



# ======================================================
# Création de l'exécutable sequentiel-1 - version 1
# ======================================================
$(BIN)/$(SQT1): $(OBJ)/$(SQT1)/affichage.o $(OBJ)/$(SQT1)/calcul-mat.o $(OBJ)/$(SQT1)/convert.o $(OBJ)/$(SQT1)/resolution.o $(OBJ)/$(SQT1)/main.o
	$(CC) $(FLAGS) $^ -o $@

$(OBJ)/$(SQT1)/main.o: $(SRC)/$(SQT1)/main.c
	$(CC) $(FLAGS) -c $< -o $@

$(OBJ)/$(SQT1)/resolution.o: $(SRC)/$(SQT1)/resolution.c
	$(CC) $(FLAGS) -c $< -o $@

$(OBJ)/$(SQT1)/affichage.o: ../$(FCC)/affichage.c
	$(CC) $(FLAGS) -c $< -o $@

$(OBJ)/$(SQT1)/calcul-mat.o: ../$(FCC)/calcul-mat.c
	$(CC) $(FLAGS) -c $< -o $@

$(OBJ)/$(SQT1)/convert.o: ../$(FCC)/convert.c
	$(CC) $(FLAGS) -c $< -o $@



.PHONY: clean

all: clean sequentiel-1

sequentiel-1: $(BIN)/$(SQT1)

exe-sequentiel-1:
	make $(SQT1)
	$(BIN)/$(SQT1) 0 1 200 1 200

clean:
	rm -f ./texte/*.data ./texte/*.txt
	rm -f $(OBJ)/$(SQT1)/*.o $(BIN)/* $(TXT)/resultats.txt