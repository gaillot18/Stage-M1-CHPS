# include <stdio.h>
# include <stdlib.h>
# include <sys/time.h>
# include <math.h>
# include <float.h>
# include <limits.h>

# include "../../Librairies/parallele-1.h"

# define pi 3.14159265358979323846

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



// Initialiser u_anc
void init_u_anc(double **u_anc){

    *u_anc = (double *)malloc(nb_pt * sizeof(double));
    
    for (int i = 0 ; i < nb_pt ; i ++){
        (*u_anc)[i] = 0.0;
    }

}



// Appliquer le schéma à un point
static inline __attribute__((always_inline)) double schema(double *f, double *u_anc, int i){

    double res = 0.5 * ((u_anc[i - 1] + u_anc[i + 1]) + h_carre * f[i]);

    return res;

}



// Calculer la norme infinie relative
static inline __attribute__((always_inline, unused)) double norme_infty_iteration(double *u, double *u_anc){

    double norme_nume = 0.0;
    double norme_deno = 0.0;
    double norme;

    # pragma omp parallel for schedule (runtime)
    for (int i = 1 ; i < nb_pt ; i ++){
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
        # pragma omp parallel for schedule(runtime)
        for (int i = 1 ; i < nb_pt - 1 ; i ++){
            u[i] = schema(f, u_anc, i);
        }

        // Test d'arrêt
        norme = norme_infty_diff(u, u_anc, nb_pt) / norme_infty(u_anc, nb_pt);

        permut = u; u = u_anc; u_anc = permut; nb_iteration ++;
        
    }

    terminaison(&permut, &u, &u_anc);

}