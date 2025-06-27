# include <stdio.h>
# include <stdlib.h>
# include <omp.h>
# include <sys/time.h>
# include <math.h>
# include <float.h>
# include <limits.h>

# include "../../Librairies/sequentiel-2.h"

# define pi 3.14159265358979323846
# define IDX(i, j) ((j) * (nb_pt) + (i))
# define idx_max (((N) - 1) * ((N) - 1))
# define max(a, b) ((a) > (b) ? (a) : (b))
# define min(a, b) ((a) < (b) ? (a) : (b))

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



// Initialiser les bords de u à 0
void init_u_bord(double *u){

    for (int i = 0 ; i < nb_pt ; i ++){
        u[IDX(i, 0)] = 0.0;
    }
    for (int i = 0 ; i < nb_pt ; i ++){
        u[IDX(i, nb_pt - 1)] = 0.0;
    }
    for (int j = 0 ; j < nb_pt ; j ++){
        u[IDX(0, j)] = 0.0;
    }
    for (int j = 0 ; j < nb_pt ; j ++){
        u[IDX(nb_pt - 1, j)] = 0.0;
    }

}



// Allouer une structure mat_Nbandes
void init_mat_Nbandes(struct mat_Nbandes *A){

    A -> N = N;
    A -> diags = (double **)malloc(N * sizeof(double *));
    for (int i = 0 ; i < N ; i ++){
        (A -> diags)[i] = (double *)malloc((idx_max - i) * sizeof(double));
    }

}



// Libérer une structure mat_Nbandes
void liberer_mat_Nbandes(struct mat_Nbandes *A){

    int N = A -> N;
    for (int i = 0 ; i < N ; i ++){
        free((A -> diags)[i]);
    }
    free(A -> diags);

}



// Connaître a_{i,j} en fonction de i et de j
static inline __attribute__((always_inline)) double valeur_a(int i, int j){

    double res;

    if (i == j){
        res = 4.0 / h_carre;
    }
    else if (i == j + 1 && j % (N - 1) != (N - 2)){
        res = -1.0 / h_carre;
    }
    else if (i == j + N - 1){
        res = -1.0 / h_carre;
    }
    else{
        res = 0.0;
    }

    return res;

}



// Obtenir la décomposition de Cholesky
void calculer_cholesky(struct mat_Nbandes *L){

    h_carre = 1.0 / pow(N, 2);
    double alpha = 4.0 / h_carre;

    for (int j = 0 ; j < idx_max ; j ++){
        for (int d = 0 ; d < N && j + d < idx_max ; d ++){

            int i = d + j;

            if (d == 0){
                (L -> diags)[0][j] = alpha;
                for (int k = max(0, j - N + d + 1) ; k < i ; k ++){
                    int d_1 = i - k;
                    (L -> diags)[0][j] -= pow((L -> diags)[d_1][k], 2);
                }
                (L -> diags)[0][j] = sqrt((L -> diags)[0][j]);
            }

            else{
                (L -> diags)[d][j] = valeur_a(i, j);
                for (int k = max(0, j - N + d + 1) ; k < j ; k ++){
                    int d_1 = i - k;
                    int d_2 = j - k;
                    (L -> diags)[d][j] -= (L -> diags)[d_1][k] * (L -> diags)[d_2][k];
                }
                (L -> diags)[d][j] /= (L -> diags)[0][j];
            }

        }
    }

}



// Résoudre Ly = f (descente)
void resoudre_cholesky_descente(struct mat_Nbandes *L, double *f, double *y){

    y[0] = f[0] / (L -> diags)[0][0];

    for (int i = 1 ; i < idx_max ; i ++){
        y[i] = f[i];
        for (int k = max(0, i - N + 1) ; k < i ; k ++){
            int d = i - k;
            y[i] -= (L -> diags)[d][k] * y[k];
        }
        y[i] /= (L -> diags)[0][i];
    }

}



// Résoudre L^{T}u = y (remontée)
void resoudre_cholesky_remontee(struct mat_Nbandes *L, double *y, double *u){
    
    u[idx_max - 1] = y[idx_max - 1] / (L -> diags)[0][idx_max - 1];

    for (int i = idx_max - 2 ; i >= 0 ; i --){
        u[i] = y[i];
        for (int k = i + 1 ; k < min(i + N, idx_max) ; k ++){
            int d = k - i;
            u[i] -= (L -> diags)[d][i] * u[k];
        }
        u[i] /= (L -> diags)[0][i];
    }

}



// Résoudre Au = f avec les conditions aux bords
void resoudre_cholesky(double *f, double *u){
    
    struct mat_Nbandes L;
    double *y, *u_int, *f_int;

    init_u_bord(u);
    u_int = (double *)malloc(idx_max * sizeof(double));
    f_int = (double *)malloc(idx_max * sizeof(double));
    extraire_interieur(u, u_int, nb_pt);
    extraire_interieur(f, f_int, nb_pt);
    y = (double *)malloc(idx_max * sizeof(double));

    init_mat_Nbandes(&L);
    calculer_cholesky(&L);

    resoudre_cholesky_descente(&L, f_int, y);
    resoudre_cholesky_remontee(&L, y, u_int);
    inserer_interieur(u_int, u, nb_pt);

    liberer_mat_Nbandes(&L);
    free(u_int);
    free(f_int);
    free(y);

}