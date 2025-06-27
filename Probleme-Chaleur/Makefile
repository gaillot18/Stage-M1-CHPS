CC = clang
CCP = mpicc

CPPFLAGS = -I/opt/homebrew/opt/llvm/include -I/opt/homebrew/opt/suite-sparse/include/suitesparse
LDFLAGS = -L/opt/homebrew/opt/llvm/lib -L/opt/homebrew/opt/suite-sparse/lib 
LDLIBS = -lcholmod

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
SQT3 = sequentiel-3

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



# ======================================================
# Création de l'exécutable parallele-1 - version 2
# ======================================================
$(BIN)/$(PRL1): $(OBJ)/$(PRL1)/affichage.o $(OBJ)/$(PRL1)/calcul-mat.o $(OBJ)/$(PRL1)/convert.o $(OBJ)/$(PRL1)/resolution.o $(OBJ)/$(PRL1)/main.o
	$(CC) $(OMP_FLAGS) $(FLAGS) $^ -o $@

$(OBJ)/$(PRL1)/main.o: $(SRC)/$(PRL1)/main.c
	$(CC) $(OMP_FLAGS) $(FLAGS) -c $< -o $@

$(OBJ)/$(PRL1)/resolution.o: $(SRC)/$(PRL1)/resolution.c
	$(CC) $(OMP_FLAGS) $(FLAGS) -c $< -o $@

$(OBJ)/$(PRL1)/affichage.o: ../$(FCC)/affichage.c
	$(CC) $(FLAGS) -c $< -o $@

$(OBJ)/$(PRL1)/calcul-mat.o: ../$(FCC)/calcul-mat.c
	$(CC) $(FLAGS) -c $< -o $@

$(OBJ)/$(PRL1)/convert.o: ../$(FCC)/convert.c
	$(CC) $(FLAGS) -c $< -o $@



# ======================================================
# Création de l'exécutable parallele-2 - version 3
# ======================================================
$(BIN)/$(PRL2): $(OBJ)/$(PRL2)/affichage.o $(OBJ)/$(PRL2)/calcul-mat.o $(OBJ)/$(PRL2)/convert.o $(OBJ)/$(PRL2)/resolution.o $(OBJ)/$(PRL2)/parallele.o $(OBJ)/$(PRL2)/main.o
	$(CCP) $(MPI_FLAGS) $(FLAGS) $^ -o $@

$(OBJ)/$(PRL2)/main.o: $(SRC)/$(PRL2)/main.c
	$(CCP) $(MPI_FLAGS) $(FLAGS) -c $< -o $@

$(OBJ)/$(PRL2)/parallele.o: $(SRC)/$(PRL2)/parallele.c
	$(CCP) $(MPI_FLAGS) $(FLAGS) -c $< -o $@

$(OBJ)/$(PRL2)/resolution.o: $(SRC)/$(PRL2)/resolution.c
	$(CCP) $(MPI_FLAGS) $(FLAGS) -c $< -o $@

$(OBJ)/$(PRL2)/affichage.o: ../$(FCC)/affichage.c
	$(CCP) $(MPI_FLAGS) $(FLAGS) -c $< -o $@

$(OBJ)/$(PRL2)/calcul-mat.o: ../$(FCC)/calcul-mat.c
	$(CCP) $(MPI_FLAGS) $(FLAGS) -c $< -o $@

$(OBJ)/$(PRL2)/convert.o: ../$(FCC)/convert.c
	$(CCP) $(MPI_FLAGS) $(FLAGS) -c $< -o $@



# ======================================================
# Création de l'exécutable sequentiel-2 - version 4
# ======================================================
$(BIN)/$(SQT2): $(OBJ)/$(SQT2)/affichage.o $(OBJ)/$(SQT2)/calcul-mat.o $(OBJ)/$(SQT2)/convert.o $(OBJ)/$(SQT2)/resolution.o $(OBJ)/$(SQT2)/main.o
	$(CC) $(FLAGS) $(CPPFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(OBJ)/$(SQT2)/main.o: $(SRC)/$(SQT2)/main.c
	$(CC) $(FLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJ)/$(SQT2)/resolution.o: $(SRC)/$(SQT2)/resolution.c
	$(CC) $(FLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJ)/$(SQT2)/affichage.o: ../$(FCC)/affichage.c
	$(CC) $(FLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJ)/$(SQT2)/calcul-mat.o: ../$(FCC)/calcul-mat.c
	$(CC) $(FLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJ)/$(SQT2)/convert.o: ../$(FCC)/convert.c
	$(CC) $(FLAGS) $(CPPFLAGS) -c $< -o $@



# ======================================================
# Création de l'exécutable sequentiel-3 - version 5
# ======================================================
$(BIN)/$(SQT3): $(OBJ)/$(SQT3)/affichage.o $(OBJ)/$(SQT3)/calcul-mat.o $(OBJ)/$(SQT3)/convert.o $(OBJ)/$(SQT3)/resolution.o $(OBJ)/$(SQT3)/main.o
	$(CC) $(FLAGS) $(CPPFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(OBJ)/$(SQT3)/main.o: $(SRC)/$(SQT3)/main.c
	$(CC) $(FLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJ)/$(SQT3)/resolution.o: $(SRC)/$(SQT3)/resolution.c
	$(CC) $(FLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJ)/$(SQT3)/affichage.o: ../$(FCC)/affichage.c
	$(CC) $(FLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJ)/$(SQT3)/calcul-mat.o: ../$(FCC)/calcul-mat.c
	$(CC) $(FLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJ)/$(SQT3)/convert.o: ../$(FCC)/convert.c
	$(CC) $(FLAGS) $(CPPFLAGS) -c $< -o $@



.PHONY: clean

all: clean sequentiel-1 parallele-1 parallele-2 sequentiel-2 sequentiel-3

sequentiel-1: $(BIN)/$(SQT1)

parallele-1: $(BIN)/$(PRL1)

parallele-2: $(BIN)/$(PRL2)

sequentiel-2: $(BIN)/$(SQT2)

sequentiel-3: $(BIN)/$(SQT3)

exe-sequentiel-1:
	make clean
	make $(SQT1)
	rm -f $(TXT)/$(SQT1)/*
	$(BIN)/$(SQT1) 0 1 200 1 160000

exe-parallele-1:
	make clean
	make $(PRL1)
	rm -f $(TXT)/$(PRL1)/*
	OMP_NUM_THREADS=4 $(BIN)/$(PRL1) 0 1 200 1 160000

exe-parallele-2:
	make clean
	make $(PRL2)
	rm -f $(TXT)/$(PRL2)/*
	mpiexec -n 4 $(BIN)/$(PRL2) 0 1 200 1 160000

exe-sequentiel-2:
	make clean
	make $(SQT2)
	rm -f $(TXT)/$(SQT2)/*
	$(BIN)/$(SQT2) 0 1 1000 1 1000

exe-sequentiel-3:
	make clean
	make $(SQT3)
	rm -f $(TXT)/$(SQT3)/*
	$(BIN)/$(SQT3) 1 1 0.02 1 0.00001

clean:
	rm -f $(OBJ)/$(SQT1)/*.o $(OBJ)/$(PRL1)/*.o $(OBJ)/$(PRL2)/*.o $(OBJ)/$(SQT2)/*.o $(OBJ)/$(SQT3)/*.o $(BIN)/* $(TXT)/resultats.txt