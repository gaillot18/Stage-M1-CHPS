# include <stdio.h>
# include <stdlib.h>
# include <sys/time.h>
# include <math.h>
# include <float.h>

# include "../../Librairies/sequentiel-2.h"

# define pi 3.14159265358979323846
# define idx_max ((N) - 1)

double h_carre;



// f dont on connait la solution exacte (1)
void f_0(double **f){

    *f = (double *)malloc(nb_pt * sizeof(double));
    for (int i = 0 ; i < nb_pt ; i ++){
        (*f)[i] = 1;
    }

}



// Solution exacte (1)
double u_e_0(double x){

    double res = 0.5 * x * (1 - x);

    return res;

}



// f dont on connait la solution exacte (2)
void f_1(double **f){

    *f = (double *)malloc(nb_pt * sizeof(double));
    double h = 1.0 / N;
    for (int i = 0 ; i < nb_pt ; i ++){
        (*f)[i] = pi * pi * sin(pi * i * h);
    }

}



// Solution exacte (2)
double u_e_1(double x){

    double res = sin(pi * x);

    return res;

}



// Calculer la solution exacte
void calculer_u_exact(double (*fonction)(double), double *u){

    double h = 1.0 / N;
    for (int i = 0 ; i < nb_pt ; i ++){
        u[i] = fonction(i * h);
    }

}



// Allouer une structure mat_2bandes
void init_mat_2bandes(struct mat_2bandes *A){

    A -> N = N;
    A -> diag = (double *)malloc(idx_max * sizeof(double));
    A -> sous_diag = (double *)malloc((idx_max - 1) * sizeof(double));

}



// Libérer une structure mat_2bandes
void liberer_mat_2bandes(struct mat_2bandes *A){

    free(A -> diag);
    free(A -> sous_diag);

}



// Obtenir la décomposition de Cholesky
void calculer_cholesky(struct mat_2bandes *L){

    h_carre = 1.0 / pow(N, 2);
    double alpha = 2.0 / h_carre;
    double beta = -1.0 / h_carre;

    (L -> diag)[0] = sqrt(alpha);
    (L -> sous_diag)[0] = beta / (L -> diag)[0];

    for (int i = 1 ; i < idx_max - 1 ; i ++){
        (L -> diag)[i] = sqrt(alpha - pow((L -> sous_diag[i - 1]), 2));
        (L -> sous_diag)[i] = beta / (L -> diag[i]);
    }

    (L -> diag)[idx_max - 1] = sqrt(alpha - pow((L -> sous_diag[idx_max - 2]), 2));

}



// Résoudre Ly = f (descente)
void resoudre_cholesky_descente(struct mat_2bandes *L, double *f, double *y){

    y[0] = f[0] / (L -> diag)[0];

    for (int i = 1 ; i < idx_max ; i ++){
        y[i] = (f[i] - (L -> sous_diag)[i - 1] * y[i - 1]) / (L -> diag)[i];
    }

}



// Résoudre L^{T}u = y (remontée)
void resoudre_cholesky_remontee(struct mat_2bandes *L, double *y, double *u){

    u[idx_max - 1] = y[idx_max - 1] / (L -> diag)[idx_max - 1];

    for (int i = idx_max - 2 ; i >= 0 ; i --){
        u[i] = (y[i] - (L -> sous_diag)[i] * u[i + 1]) / (L -> diag)[i];
    }

}



// Résoudre Au = f avec les conditions aux bords
void resoudre_cholesky(double *f, double *u){
    
    struct mat_2bandes L;
    double *y = (double *)malloc(idx_max * sizeof(double));

    u[0] = 0; u[nb_pt - 1] = 0;

    init_mat_2bandes(&L);
    calculer_cholesky(&L);

    resoudre_cholesky_descente(&L, &(f[1]), y); // Laisser f[0] pour le bord
    resoudre_cholesky_remontee(&L, y, &(u[1])); // Laisser u[0] pour le bord

    liberer_mat_2bandes(&L);
    free(y);

}