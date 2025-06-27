# include <stdio.h>
# include <stdlib.h>
# include <sys/time.h>
# include <math.h>
# include <float.h>
# include <limits.h>
# include <cholmod.h>

# include "../../Librairies/sequentiel-3.h"

# define pi 3.14159265358979323846
# define IDX(i, j) ((j) * (nb_pt) + (i))
# define idx_max (((N) - 1) * ((N) - 1))
# define nb_elements (((N) - 3) * (5 * (N) + 1) + 12)

double h_carre;



// f dont on connait la solution exacte
void f_1(double **f){

    *f = (double *)malloc(nb_pt * nb_pt * sizeof(double));
    double h = 1.0 / N;
    for (int j = 0 ; j < nb_pt ; j ++){
        for (int i = 0 ; i < nb_pt ; i ++){
            (*f)[IDX(i, j)] = sin(2 * pi * i * h) * sin(2 * pi * j * h);
        }
    }

}



// Solution exacte
double u_e_1(double x, double y){

    double res = 1.0 / (8 * pow(pi, 2)) * sin(2 * pi * x) * sin(2 * pi * y);

    return res;

}



// Calculer la solution exacte
void calculer_u_exact(double (*fonction)(double, double), double *u){

    double h = 1.0 / N;
    for (int j = 0 ; j < nb_pt ; j ++){
        for (int i = 0 ; i < nb_pt ; i ++){
            u[IDX(i, j)] = fonction(i * h, j * h);
        }
    }

}



// Connaitre le type de bord (voir figure du rapport)
static inline __attribute__((always_inline)) int connaitre_bord(int x, int y){

    int res;

    if (x == 0){
        if (y == 0){res = 1;}
        else if (y == N - 2){res = 2;}
        else{res = -1;}
    }
    else if (y == 0){
        if (x == N - 2){res = 4;}
        else{res = -4;}
    }
    else if (x == N - 2){
        if (y == N - 2){res = 3;}
        else{res = -3;}
    }
    else if (y == N - 2){res = - 2;}
    else{res = 0;}

    return res;

}



// Construire les tableaux lignes, valeurs et offsets (format CSC) pour avoir la matrice creuse de A
void construire_matrice_creuse(int **lignes, double **valeurs, int **offsets){

    h_carre = 1.0 / pow(N, 2);

    double val_coin_bas_gauche[3] = {4 / h_carre, -1 / h_carre, -1 / h_carre}; // 1
    int lig_coin_bas_gauche[3] = {0, 1, N - 1};

    double val_coin_haut_droite[3] = {-1 / h_carre, -1 / h_carre, 4 / h_carre}; // 3
    int lig_coin_haut_droite[3] = {-N + 1, -1, 0};

    double val_coin_bas_droite[3] = {-1 / h_carre, 4 / h_carre, -1 / h_carre}; // 4
    int lig_coin_bas_droite[3] = {-1, 0, N - 1};

    double val_coin_haut_gauche[3] = {-1 / h_carre, 4 / h_carre, -1 / h_carre}; // 2
    int lig_coin_haut_gauche[3] = {-N + 1, 0, 1};

    double val_bord_bas[4] = {-1 / h_carre, 4 / h_carre, -1 / h_carre, -1 / h_carre}; // -4
    int lig_bord_bas[4] = {-1, 0, 1, N - 1};

    double val_bord_gauche[4] = {-1 / h_carre, 4 / h_carre, -1 / h_carre, -1 / h_carre}; // -1
    int lig_bord_gauche[4] = {-N + 1, 0, 1, N - 1};

    double val_bord_haut[4] = {-1 / h_carre, -1 / h_carre, 4 / h_carre, -1 / h_carre}; // -2
    int lig_bord_haut[4] = {-N + 1, -1, 0, 1};

    double val_bord_droite[4] = {-1 / h_carre, -1 / h_carre, 4 / h_carre, -1 / h_carre}; // -3
    int lig_bord_droite[4] = {-N + 1, -1, 0, N - 1};

    double val_interieur[5] = {-1 / h_carre, -1 / h_carre, 4 / h_carre, -1 / h_carre, -1 / h_carre}; // 0
    int lig_interieur[5] = {-N + 1, -1, 0, 1, N - 1};

    *lignes = (int *)malloc(nb_elements * sizeof(int));
    *valeurs = (double *)malloc(nb_elements * sizeof(double));
    *offsets = (int *)malloc((idx_max + 1) * sizeof(int));
    int *source_ligne; double *source_valeur;
    int idx = 0; int offset;
    (*offsets)[0] = 0;

    for (int j = 0 ; j < idx_max ; j ++){
        int x = j % (N - 1);
        int y = j / (N - 1);
        int bord = connaitre_bord(x, y);
        if (bord == 0){
            source_ligne = lig_interieur;
            source_valeur = val_interieur;
            offset = 5;
        }
        else if (bord == -1){
            source_ligne = lig_bord_gauche;
            source_valeur = val_bord_gauche;
            offset = 4;
        }
        else if (bord == -2){
            source_ligne = lig_bord_haut;
            source_valeur = val_bord_haut;
            offset = 4;
        }
        else if (bord == -3){
            source_ligne = lig_bord_droite;
            source_valeur = val_bord_droite;
            offset = 4;
        }
        else if (bord == -4){
            source_ligne = lig_bord_bas;
            source_valeur = val_bord_bas;
            offset = 4;
        }
        else if (bord == 1){
            source_ligne = lig_coin_bas_gauche;
            source_valeur = val_coin_bas_gauche;
            offset = 3;
        }
        else if (bord == 2){
            source_ligne = lig_coin_haut_gauche;
            source_valeur = val_coin_haut_gauche;
            offset = 3;
        }
        else if (bord == 3){
            source_ligne = lig_coin_haut_droite;
            source_valeur = val_coin_haut_droite;
            offset = 3;
        }
        else{
            source_ligne = lig_coin_bas_droite;
            source_valeur = val_coin_bas_droite;
            offset = 3;
        }

        for (int i = 0 ; i < offset ; i ++){
            (*lignes)[idx] = j + source_ligne[i];
            (*valeurs)[idx] = source_valeur[i];
            idx ++;
        }
        (*offsets)[j + 1] = (*offsets)[j] + offset;
    }

}



// Allouer et initialiser la matrice creuse A
cholmod_sparse *init_matrice_creuse(int *offsets, int *lignes, double *valeurs){

    cholmod_sparse *A = cholmod_allocate_sparse(idx_max, idx_max, nb_elements, 1, 1, 0, CHOLMOD_REAL, &c);

    memcpy(A -> p, offsets, (idx_max + 1) * sizeof(int));
    memcpy(A -> i, lignes, nb_elements * sizeof(int));
    memcpy(A -> x, valeurs, nb_elements * sizeof(double));
    A -> nzmax = nb_elements;
    A -> stype = 1;
    A -> nrow = idx_max;
    A -> ncol = idx_max;

    free(lignes);
    free(valeurs);
    free(offsets);

    return A;

}



// Fonction principale
void resoudre(cholmod_sparse *A, double *f, double *u){

    h_carre = 1.0 / pow(N, 2);
    double *f_int = (double *)malloc(idx_max * sizeof(double));
    double *u_int = (double *)malloc(idx_max * sizeof(double));
    for (int i = 0 ; i < nb_pt * nb_pt ; i ++){
        u[i] = 0.0;
    }

    extraire_interieur(f, f_int, nb_pt);
    extraire_interieur(u, u_int, nb_pt);

    cholmod_dense *f_dense = cholmod_allocate_dense(A -> nrow, 1, A -> nrow, CHOLMOD_REAL, &c);
    memcpy(f_dense -> x, f_int, A -> nrow * sizeof(double));

    cholmod_factor *L = cholmod_analyze(A, &c);
    cholmod_factorize(A, L, &c);

    cholmod_dense *u_dense = cholmod_solve(CHOLMOD_A, L, f_dense, &c);

    memcpy(u_int, u_dense -> x, A -> nrow * sizeof(double));

    inserer_interieur(u_int, u, nb_pt);

    cholmod_free_factor(&L, &c);
    cholmod_free_dense(&f_dense, &c);
    cholmod_free_dense(&u_dense, &c);
    free(f_int);
    free(u_int);

}