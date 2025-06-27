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

.SUFFIXES: .o .c



# ======================================================
# Création de l'exécutable base - version 0
# ======================================================
$(BIN)/$(BS): $(OBJ)/$(BS)/affichage.o $(OBJ)/$(BS)/calcul-mat.o $(OBJ)/$(BS)/convert.o $(OBJ)/$(BS)/resolution.o $(OBJ)/$(BS)/main.o
	$(CC) $(FLAGS) $^ -o $@

$(OBJ)/$(BS)/main.o: $(SRC)/$(BS)/main.c
	$(CC) $(FLAGS) -c $< -o $@

$(OBJ)/$(BS)/resolution.o: $(SRC)/$(BS)/resolution.c
	$(CC) $(FLAGS) -c $< -o $@

$(OBJ)/$(BS)/affichage.o: ../$(FCC)/affichage.c
	$(CC) $(FLAGS) -c $< -o $@

$(OBJ)/$(BS)/calcul-mat.o: ../$(FCC)/calcul-mat.c
	$(CC) $(FLAGS) -c $< -o $@

$(OBJ)/$(BS)/convert.o: ../$(FCC)/convert.c
	$(CC) $(FLAGS) -c $< -o $@



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
$(BIN)/$(SQT2): $(OBJ)/$(SQT2)/affichage.o $(OBJ)/$(SQT2)/calcul-mat.o $(OBJ)/$(SQT2)/convert.o $(OBJ)/$(SQT2)/resolution.o $(OBJ)/$(SQT2)/affichage-struct.o $(OBJ)/$(SQT2)/main.o
	$(CC) $(FLAGS) $^ -o $@

$(OBJ)/$(SQT2)/main.o: $(SRC)/$(SQT2)/main.c
	$(CC) $(FLAGS) -c $< -o $@

$(OBJ)/$(SQT2)/resolution.o: $(SRC)/$(SQT2)/resolution.c
	$(CC) $(FLAGS) -c $< -o $@

$(OBJ)/$(SQT2)/affichage-struct.o: $(SRC)/$(SQT2)/affichage-struct.c
	$(CC) $(FLAGS) -c $< -o $@

$(OBJ)/$(SQT2)/affichage.o: ../$(FCC)/affichage.c
	$(CC) $(FLAGS) -c $< -o $@

$(OBJ)/$(SQT2)/calcul-mat.o: ../$(FCC)/calcul-mat.c
	$(CC) $(FLAGS) -c $< -o $@

$(OBJ)/$(SQT2)/convert.o: ../$(FCC)/convert.c
	$(CC) $(FLAGS) -c $< -o $@



.PHONY: clean

all: clean base sequentiel-1 parallele-1 parallele-2 sequentiel-2

base : $(BIN)/$(BS)

sequentiel-1: $(BIN)/$(SQT1)

parallele-1: $(BIN)/$(PRL1)

parallele-2: $(BIN)/$(PRL2)

sequentiel-2: $(BIN)/$(SQT2)

exe-base:
	make $(BS)
	$(BIN)/$(BS) 100

exe-sequentiel-1:
	make $(SQT1)
	$(BIN)/$(SQT1) 100

exe-parallele-1:
	make $(PRL1)
	OMP_NUM_THREADS=4 $(BIN)/$(PRL1) 100

exe-parallele-2:
	make $(PRL2)
	mpiexec -n 4 $(BIN)/$(PRL2) 100

exe-sequentiel-2:
	make $(SQT2)
	$(BIN)/$(SQT2) 100

clean:
	rm -f $(OBJ)/$(BS)/*.o $(OBJ)/$(SQT1)/*.o $(OBJ)/$(PRL1)/*.o $(OBJ)/$(PRL2)/*.o $(OBJ)/$(SQT2)/*.o $(BIN)/* $(TXT)/resultats.txt