# include <stdio.h>
# include <stdlib.h>
# include <omp.h>
# include <sys/time.h>
# include <math.h>
# include <float.h>
# include <limits.h>

# include "../../Librairies/base.h"

# define pi 3.14159265358979323846
# define idx_max (((N) - 1))
# define IDX(i, j) ((j) * (nb_pt) + (i))
# define IDX_M(i, j) ((i) * ((N) - 1) + (j))

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



// Construire la matrice A
double *construire_matrice(){

    h_carre = 1.0 / pow(N, 2);
    double alpha = 2.0 / h_carre;
    double beta = -1.0 / h_carre;
    double *A = (double *)malloc(idx_max * idx_max * sizeof(double));
    for (int i = 0 ; i < idx_max * idx_max ; i ++){
        A[i] = 0.0;
    }

    A[IDX_M(0, 0)] = alpha; A[IDX_M(1, 0)] = beta;
    A[IDX_M(idx_max - 1, idx_max - 1)] = alpha; A[IDX_M(idx_max - 2, idx_max - 1)] = beta;

    for (int j = 1 ; j < idx_max - 1 ; j ++){

        A[IDX_M(j - 1, j)] = beta;
        A[IDX_M(j, j)] = alpha;
        A[IDX_M(j + 1, j)] = beta;

    }

    return A;

}



// Fonction principale
void resoudre_gauss(double *A, double *f, double *u){

    double *f_int = (double *)malloc(idx_max * sizeof(double));
    double *u_int = (double *)malloc(idx_max * sizeof(double));
    for (int i = 0 ; i < idx_max ; i ++){
        f_int[i] = f[i + 1];
        u_int[i] = u[i + 1];
    }

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

    u[0] = 0.0; u[nb_pt - 1] = 0.0;
    for (int i = 1 ; i < nb_pt - 1 ; i ++){
        u[i] = u_int[i - 1];
    }

    free(f_int);
    free(u_int);

}