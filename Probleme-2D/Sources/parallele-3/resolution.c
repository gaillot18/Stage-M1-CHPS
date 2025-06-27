# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <sys/time.h>
# include <math.h>
# include <mpi.h>
# include <float.h>
# include <limits.h>

# include "../../Librairies/parallele-3.h"

# define pi 3.14159265358979323846
# define IDX(i, j) ((j) * (nb_pt_div_i + 2) + (i))

double h_carre;



// f dont on connait la solution exacte
void f_1(double **f){

    *f = (double *)malloc((nb_pt_div_i + 2) * (nb_pt_div_j + 2) * sizeof(double));
    double h = 1.0 / N;
    int j_reel = j_debut;
    int i_reel;

    for (int i = 0 ; i < (nb_pt_div_i + 2) * (nb_pt_div_j + 2) ; i ++){ // A améliorer
        (*f)[i] = 0.0;
    }
    
    for (int j = 1 ; j < nb_pt_div_j + 1 ; j ++){
        i_reel = i_debut;
        for (int i = 1 ; i < nb_pt_div_i + 1 ; i ++){
            (*f)[IDX(i, j)] = sin(2 * pi * i_reel * h) * sin(2 * pi * j_reel * h);
            i_reel ++;
        }
        j_reel ++;
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
            u[j * nb_pt + i] = fonction(i * h, j * h);
        }
    }

}



// Initialiser u_div_anc
void init_u_anc(double **u_div_anc){

    *u_div_anc = (double *)malloc((nb_pt_div_i + 2) * (nb_pt_div_j + 2) * sizeof(double));

    for (int i = 0 ; i < (nb_pt_div_i + 2) * (nb_pt_div_j + 2) ; i ++){
        (*u_div_anc)[i] = 0.0;
    }

}



// Appliquer le schéma à un point
static inline __attribute__((always_inline)) double schema(double *f_div, double *u_div_anc, int i, int j){

    double res = 0.25 * (
    u_div_anc[IDX(i - 1, j)]
    + u_div_anc[IDX(i, j - 1)]
    + u_div_anc[IDX(i + 1, j)]
    + u_div_anc[IDX(i, j + 1)]
    + h_carre * f_div[IDX(i, j)]);

    return res;

}



// Appliquer le schéma sur les bords
static inline __attribute__((always_inline)) void calculer_u_jacobi_bords(double *f_div, double *u_div, double *u_div_anc){

    // Bords sauf coins
    if (j_debut != 0){ 
        int j = 1;
        for (int i = 2 ; i < nb_pt_div_i ; i ++){
            u_div[IDX(i, j)] = schema(f_div, u_div_anc, i, j);
        }
    }

    if (j_fin != nb_pt - 1){
        int j = nb_pt_div_j;
        for (int i = 2 ; i < nb_pt_div_i ; i ++){
            u_div[IDX(i, j)] = schema(f_div, u_div_anc, i, j);
        }
    }

    if (i_debut != 0){
        int i = 1;
        for (int j = 2 ; j < nb_pt_div_j ; j ++){
            u_div[IDX(i, j)] = schema(f_div, u_div_anc, i, j);
        }
    }
    
    if (i_fin != nb_pt - 1){
        int i = nb_pt_div_i;
        for (int j = 2 ; j < nb_pt_div_j ; j ++){
            u_div[IDX(i, j)] = schema(f_div, u_div_anc, i, j);
        }
    }

    // Coins
    if (i_debut != 0 && j_debut != 0){
        int i = 1; int j = 1;
        u_div[IDX(i, j)] = schema(f_div, u_div_anc, i, j);
    }

    if (i_fin != nb_pt - 1 && j_debut != 0){
        int i = nb_pt_div_i; int j = 1;
        u_div[IDX(i, j)] = schema(f_div, u_div_anc, i, j);
    }

    if (i_fin != nb_pt - 1 && j_fin != nb_pt - 1){
        int i = nb_pt_div_i; int j = nb_pt_div_j;
        u_div[IDX(i, j)] = schema(f_div, u_div_anc, i, j);
    }

    if (i_debut != 0 && j_fin != nb_pt - 1){
        int i = 1; int j = nb_pt_div_j;
        u_div[IDX(i, j)] = schema(f_div, u_div_anc, i, j);
    }

}



// Calculer la norme infinie relative
static inline __attribute__((always_inline)) double norme_infty_iteration(double *u_div, double *u_div_anc){

    double norme_nume_div = 0.0;
    double norme_deno_div = 0.0;
    double norme_nume;
    double norme_deno;
    double norme;

    for (int j = 1 ; j < nb_pt_div_j + 1 ; j ++){
        for (int i = 1 ; i < nb_pt_div_i + 1 ; i ++){
            double diff = fabs(u_div[IDX(i, j)] - u_div_anc[IDX(i, j)]);
            if (diff > norme_nume_div){
                norme_nume_div = diff;
            }
            if (fabs(u_div_anc[IDX(i, j)]) > norme_deno_div){
                norme_deno_div = fabs(u_div_anc[IDX(i, j)]);
            }
        }
    }

    MPI_Allreduce(&norme_nume_div, &norme_nume, 1, MPI_DOUBLE, MPI_MAX, comm_2D);
    MPI_Allreduce(&norme_deno_div, &norme_deno, 1, MPI_DOUBLE, MPI_MAX, comm_2D);
    norme = norme_nume / norme_deno;

    return norme;

}



// Terminer
void terminaison(double **permut, double **u_div, double **u_div_anc){

    if (nb_iteration % 2 != 0){
        *permut = *u_div; *u_div = *u_div_anc; *u_div_anc = *permut;
    }

    free(*u_div_anc);

}



// Fonction principale
void calculer_u_jacobi(double *f_div, double *u_div){

    nb_iteration = 0;
    h_carre = 1.0 / pow(N, 2);
    int nb_iteration_max = INT_MAX;
    double norme = DBL_MAX;
    int i_boucle_debut; int j_boucle_debut;
    int i_boucle_fin; int j_boucle_fin;
    double *u_div_anc; double *permut;

    // Vecteur de départ
    init_u_anc(&u_div_anc);
    for (int i = 0 ; i < (nb_pt_div_i + 2) * (nb_pt_div_j + 2) ; i ++){
        u_div[i] = 0.0;
    }

    // Bornes des boucles
    infos_bornes_boucles(&i_boucle_debut, &j_boucle_debut, &i_boucle_fin, &j_boucle_fin);

    // Itérations
    for (int iteration = 0 ; iteration < nb_iteration_max && norme > 1e-10 ; iteration ++){

        // Communication
        echanger_halos(u_div_anc);

        // Schéma
        for (int j = 2 ; j < nb_pt_div_j ; j ++){
            for (int i = 2 ; i < nb_pt_div_i ; i ++){
                    u_div[IDX(i, j)] = schema(f_div, u_div_anc, i, j);
            }
        }

        // Test de communication
        test_fin_echange_halos();

        // Schéma sur les bords
        calculer_u_jacobi_bords(f_div, u_div, u_div_anc);

        // Test d'arrêt
        norme = norme_infty_iteration(u_div, u_div_anc);

        permut = u_div; u_div = u_div_anc; u_div_anc = permut; nb_iteration ++;

    }

    terminaison(&permut, &u_div, &u_div_anc);

}