# include <stdio.h>
# include <stdlib.h>
# include <omp.h>
# include <sys/time.h>
# include <math.h>
# include <float.h>
# include <limits.h>

# include "../../Librairies/base.h"

# define pi 3.14159265358979323846
# define idx_max (((N) - 1) * ((N) - 1))
# define IDX(i, j) ((j) * (nb_pt) + (i))
# define IDX_M(i, j) ((i) * (((N) - 1) * ((N) - 1)) + (j))

double h_carre;


// f dont on connait la solution exacte
void f_1(double **f){

    *f = (double *)malloc(nb_pt * nb_pt * sizeof(double));
    double h = 1.0 / N;
    for (int i = 0 ; i < nb_pt ; i ++){
        for (int j = 0 ; j < nb_pt ; j ++){
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
    for (int i = 0 ; i < nb_pt ; i ++){
        for (int j = 0 ; j < nb_pt ; j ++){
            u[IDX(i, j)] = fonction(i * h, j * h);
        }
    }

}



// Initialiser u_anc
void init_u_anc(double **u_anc){

    *u_anc = (double *)malloc(nb_pt * nb_pt * sizeof(double));

    for (int i = 0 ; i < nb_pt * nb_pt ; i ++){
        (*u_anc)[i] = 0.0;
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



// Construire la matrice A
double *construire_matrice(){

    h_carre = 1.0 / pow(N, 2);
    double *A = (double *)malloc(idx_max * idx_max * sizeof(double));
    for (int i = 0 ; i < idx_max * idx_max ; i ++){
        A[i] = 0.0;
    }

    double alpha = 4.0 / h_carre;
    double beta = -1.0 / h_carre;

    double val_coin_bas_gauche[3] = {alpha, beta, beta}; // 1
    int lig_coin_bas_gauche[3] = {0, 1, N - 1};

    double val_coin_haut_droite[3] = {beta, beta, alpha}; // 3
    int lig_coin_haut_droite[3] = {-N + 1, -1, 0};

    double val_coin_bas_droite[3] = {beta, alpha, beta}; // 4
    int lig_coin_bas_droite[3] = {-1, 0, N - 1};

    double val_coin_haut_gauche[3] = {beta, alpha, beta}; // 2
    int lig_coin_haut_gauche[3] = {-N + 1, 0, 1};

    double val_bord_bas[4] = {beta, alpha, beta, beta}; // -4
    int lig_bord_bas[4] = {-1, 0, 1, N - 1};

    double val_bord_gauche[4] = {beta, alpha, beta, beta}; // -1
    int lig_bord_gauche[4] = {-N + 1, 0, 1, N - 1};

    double val_bord_haut[4] = {beta, beta, alpha, beta}; // -2
    int lig_bord_haut[4] = {-N + 1, -1, 0, 1};

    double val_bord_droite[4] = {beta, beta, alpha, beta}; // -3
    int lig_bord_droite[4] = {-N + 1, -1, 0, N - 1};

    double val_interieur[5] = {beta, beta, alpha, beta, beta}; // 0
    int lig_interieur[5] = {-N + 1, -1, 0, 1, N - 1};

    int *source_ligne; double *source_valeur;
    int offset;

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
            int ligne = j + source_ligne[i];
            A[IDX_M(ligne, j)] = source_valeur[i];
        }

    }

    return A;

}



// Fonction principale
void resoudre_gauss(double *A, double *f, double *u){

    double *f_int = (double *)malloc(idx_max * sizeof(double));
    double *u_int = (double *)malloc(idx_max * sizeof(double));
    extraire_interieur(f, f_int, nb_pt);
    extraire_interieur(u, u_int, nb_pt);

    for (int j = 0 ; j < idx_max ; j ++){
        for (int i = 0 ; i < idx_max ; i ++){
            if (i != j){
                double facteur = A[IDX_M(i, j)] / A[IDX_M(j, j)];
                for (int k = 0 ; k < idx_max ; k ++){
                    A[IDX_M(i, k)] -= facteur * A[IDX_M(j, k)];
                }
                f_int[i] -= facteur * f_int[j];
            }
        }
    }

    for (int i = idx_max - 1 ; i >= 0 ; i --){
        u_int[i] = f_int[i];
        for (int j = i + 1 ; j < idx_max ; j ++){
            u_int[i] -= A[IDX_M(i, j)] * u_int[j];
        }
        u_int[i] /= A[IDX_M(i, i)];
    }

    inserer_interieur(u_int, u, nb_pt);

    free(f_int);
    free(u_int);

}