# include <stdio.h>
# include <stdlib.h>
# include <omp.h>
# include <sys/time.h>
# include <math.h>
# include <float.h>
# include <limits.h>

# include "../../Librairies/sequentiel-1.h"

# define pi 3.14159265358979323846
# define IDX(i, j) ((j) * (nb_pt) + (i))

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



// Appliquer le schéma à un point
static inline __attribute__((always_inline)) double schema(double *f, double *u_anc, int i, int j){

    double res = 0.25 * (
    u_anc[IDX(i - 1, j)]
    + u_anc[IDX(i, j - 1)]
    + u_anc[IDX(i + 1, j)]
    + u_anc[IDX(i, j + 1)]
    + h_carre * f[IDX(i, j)]);

    return res;

}



// Calculer la norme infinie relative
static inline __attribute__((always_inline)) double norme_infty_iteration(double *u, double *u_anc){

    double norme_nume = 0.0;
    double norme_deno = 0.0;
    double norme;

    for (int i = 0 ; i < nb_pt * nb_pt ; i ++){
        double diff = fabs(u[i] - u_anc[i]);
        if (diff > norme_nume){
            norme_nume = diff;
        }
        if (fabs(u_anc[i]) > norme_deno){
            norme_deno = fabs(u_anc[i]);
        }
    }

    norme = norme_nume / norme_deno;

    return norme;

}



// Terminer
void terminaison(double **permut, double **u, double **u_anc){

    if (nb_iteration % 2 != 0){
        *permut = *u; *u = *u_anc; *u_anc = *permut;
    }

    free(*u_anc);

}



// Fonction principale
void calculer_u_jacobi(double *f, double *u){

    nb_iteration = 0;
    h_carre = 1.0 / pow(N, 2);
    int nb_iteration_max = INT_MAX;
    double norme = DBL_MAX;
    double *u_anc; double *permut;

    // Vecteur de départ
    init_u_anc(&u_anc);

    // Itérations
    for (int iteration = 0 ; iteration < nb_iteration_max && norme > 1e-10 ; iteration ++){

        // Schéma
        for (int j = 1 ; j < nb_pt - 1 ; j ++){
            for (int i = 1 ; i < nb_pt - 1 ; i ++){
                u[IDX(i, j)] = schema(f, u_anc, i, j);
            }
        }

        // Test d'arrêt
        norme = norme_infty_iteration(u, u_anc);
        
        permut = u; u = u_anc; u_anc = permut; nb_iteration ++;

    }

    terminaison(&permut, &u, &u_anc);

}