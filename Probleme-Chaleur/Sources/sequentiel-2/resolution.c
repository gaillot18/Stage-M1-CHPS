# include <stdio.h>
# include <stdlib.h>
# include <sys/time.h>
# include <math.h>
# include <float.h>
# include <limits.h>
# include <cholmod.h>

# include "../../Librairies/sequentiel-2.h"

# define pi 3.14159265358979323846
# define IDX(i, j) ((j) * (nb_pt) + (i))
# define idx_max (((N) - 1) * ((N) - 1))
# define nb_elements (((N) - 3) * (5 * (N) + 1) + 12)



// f dont on connait la solution exacte
static inline __attribute__((always_inline)) double f_source(double x, double y, double t){

    double res = (-lambda + 2 * a * pow(pi, 2)) * sin(pi * x) * sin(pi * y) * exp(-lambda * t);

    return res;

}



// Calculer le second membre b
static inline __attribute__((always_inline)) void calculer_b(int k, double *u, double *b){

    for (int i = 0 ; i < idx_max ; i ++){
        double x = (i % (N - 1) + 1) * h;
        double y = (i / (N - 1) + 1) * h;
        double t = (double)k / N_t;
        b[i] = u[i] + h_t * f_source(x, y, t + h_t);
    }

}



// Solution exacte
double u_e(double x, double y, double t){

    double res = sin(pi * x) * sin(pi * y) * exp(-lambda * t);

    return res;

}



// Calculer la solution exacte
void calculer_u_exact(double (*fonction)(double, double, double), double *u, int k){

    for (int j = 0 ; j < nb_pt ; j ++){
        for (int i = 0 ; i < nb_pt ; i ++){
            u[IDX(i, j)] = fonction(i * h, j * h, k * h_t);
        }
    }

}



// Calculer u_0
double u_0(double x, double y){

    double res = sin(pi * x) * sin(pi * y);

    return res;

}



// Allouer et initialiser u_0
void init_u_0(double (*fonction)(double, double), double *u){

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



// Écrire dans un fichier
static inline __attribute__((always_inline, unused)) void ecrire_double_iteration(double *u){

    fwrite(u, sizeof(double), nb_pt * nb_pt, descripteur);

}



// Fonction principale (calcul uniquement de la solution approchée)
void resoudre(cholmod_sparse *A, double *u){

    double *b_int = (double *)malloc(idx_max * sizeof(double));
    double *u_int = (double *)malloc(idx_max * sizeof(double));

    init_u_0(u_0, u);

    cholmod_dense *b_dense = cholmod_allocate_dense(A -> nrow, 1, A -> nrow, CHOLMOD_REAL, &c);
    cholmod_dense *u_dense;
    cholmod_factor *L = cholmod_analyze(A, &c);
    cholmod_factorize(A, L, &c);

    for (int k = 1 ; k <= N_t ; k ++){

        # ifdef ECRITURE
        ecrire_double_iteration(u);
        # endif

        extraire_interieur(u, u_int, nb_pt);
        calculer_b(k - 1, u_int, b_int);

        memcpy(b_dense -> x, b_int, A -> nrow * sizeof(double));

        u_dense = cholmod_solve(CHOLMOD_A, L, b_dense, &c);

        memcpy(u_int, u_dense -> x, A -> nrow * sizeof(double));

        inserer_interieur(u_int, u, nb_pt);

    }

    # ifdef ECRITURE
    ecrire_double_iteration(u);
    # endif

    cholmod_free_factor(&L, &c);
    cholmod_free_dense(&b_dense, &c);
    cholmod_free_dense(&u_dense, &c);
    free(b_int);
    free(u_int);


}



// Fonction principale (calcul de la solution exacte et de la solution approchée pour avoir l'erreur à chaque itération)
double resoudre_calculer_u_exact(cholmod_sparse *A, double *u){

    double *u_exact = (double *)malloc(nb_pt * nb_pt * sizeof(double));
    double erreur_infty_k; double erreur_infty = 0.0;
    double *b_int = (double *)malloc(idx_max * sizeof(double));
    double *u_int = (double *)malloc(idx_max * sizeof(double));

    init_u_0(u_0, u);

    cholmod_dense *b_dense = cholmod_allocate_dense(A -> nrow, 1, A -> nrow, CHOLMOD_REAL, &c);
    cholmod_dense *u_dense;
    cholmod_factor *L = cholmod_analyze(A, &c);
    cholmod_factorize(A, L, &c);

    for (int k = 1 ; k <= N_t ; k ++){

        extraire_interieur(u, u_int, nb_pt);
        calculer_b(k - 1, u_int, b_int);

        memcpy(b_dense -> x, b_int, A -> nrow * sizeof(double));

        u_dense = cholmod_solve(CHOLMOD_A, L, b_dense, &c);

        memcpy(u_int, u_dense -> x, A -> nrow * sizeof(double));

        inserer_interieur(u_int, u, nb_pt);

        calculer_u_exact(u_e, u_exact, k);
        erreur_infty_k = norme_infty_diff(u_exact, u, nb_pt * nb_pt);
        if (erreur_infty_k > erreur_infty){
            erreur_infty = erreur_infty_k;
        }

    }

    cholmod_free_factor(&L, &c);
    cholmod_free_dense(&b_dense, &c);
    cholmod_free_dense(&u_dense, &c);
    free(b_int);
    free(u_int);
    free(u_exact);

    return erreur_infty;

}