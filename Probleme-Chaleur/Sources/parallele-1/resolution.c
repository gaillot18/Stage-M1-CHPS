# include <stdio.h>
# include <stdlib.h>
# include <sys/time.h>
# include <math.h>
# include <float.h>
# include <limits.h>

# include "../../Librairies/parallele-1.h"

# define pi 3.14159265358979323846
# define IDX(i, j) ((j) * (nb_pt) + (i))



// f dont on connait la solution exacte
static inline __attribute__((always_inline)) double f_source(double x, double y, double t){

    double res = (-lambda + 2 * a * pow(pi, 2)) * sin(pi * x) * sin(pi * y) * exp(-lambda * t);

    return res;

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
void init_u_0(double (*fonction)(double, double), double **u_anc){

    *u_anc = (double *)malloc(nb_pt * nb_pt * sizeof(double));

    for (int j = 0 ; j < nb_pt ; j ++){
        for (int i = 0 ; i < nb_pt ; i ++){
            (*u_anc)[IDX(i, j)] = fonction(i * h, j * h);
        }
    }

}



// Appliquer le schéma à un point
static inline __attribute__((always_inline)) double schema(double f, double *u_anc, int i, int j, int k){

    double res = alpha * u_anc[IDX(i, j)]
    + beta * (u_anc[IDX(i - 1, j)] + u_anc[IDX(i, j - 1)] + u_anc[IDX(i + 1, j)] + u_anc[IDX(i, j + 1)])
    + h_t * f;

    return res;

}



// Écrire dans un fichier
static inline __attribute__((always_inline, unused)) void ecrire_double_iteration(double *u){

    fwrite(u, sizeof(double), nb_pt * nb_pt, descripteur);

}



// Terminer
void terminaison(double **permut, double **u, double **u_anc){

    if (N_t % 2 != 0){
        *permut = *u; *u = *u_anc; *u_anc = *permut;
    }

    # pragma omp single
    free(*u_anc);

}



// Fonction principale (calcul uniquement de la solution approchée)
void calculer_u(double *u){

    double *u_anc; double *permut;

    // Vecteur de départ
    init_u_0(u_0, &u_anc);
    for (int i = 0 ; i < nb_pt * nb_pt ; i ++){
        u[i] = 0.0;
    }

    # pragma omp parallel firstprivate(u, u_anc, permut)
    {

        for (int k = 1 ; k <= N_t ; k ++){

            # ifdef ECRITURE
            # pragma omp single
            ecrire_double_iteration(u_anc);
            # endif

            // Schéma
            # pragma omp for schedule(static)
            for (int j = 1 ; j < nb_pt - 1 ; j ++){
                for (int i = 1 ; i < nb_pt - 1 ; i ++){
                    double f = f_source(i * h, j * h, k * h_t);
                    u[IDX(i, j)] = schema(f, u_anc, i, j, k);
                }
            }

            permut = u; u = u_anc; u_anc = permut;

        }

        # ifdef ECRITURE
        ecrire_double_iteration(u);
        # endif

        terminaison(&permut, &u, &u_anc);

    }

}



// Fonction principale (calcul de la solution exacte et de la solution approchée pour avoir l'erreur à chaque itération)
double calculer_u_u_exact(double *u){

    double *u_exact = (double *)malloc(nb_pt * nb_pt * sizeof(double));
    double *u_anc; double *permut;
    double erreur_infty_k; double erreur_infty = 0.0;

    // Vecteur de départ
    init_u_0(u_0, &u_anc);
    for (int i = 0 ; i < nb_pt * nb_pt ; i ++){
        u[i] = 0.0;
    }

    # pragma omp parallel firstprivate(u, u_anc, permut)
    {

        for (int k = 1 ; k <= N_t ; k ++){

            // Schéma
            # pragma omp for schedule(static)
            for (int j = 1 ; j < nb_pt - 1 ; j ++){
                for (int i = 1 ; i < nb_pt - 1 ; i ++){
                    double f = f_source(i * h, j * h, k * h_t);
                    u[IDX(i, j)] = schema(f, u_anc, i, j, k);
                }
            }

            // Calcul de la solution exacte
            # pragma omp single
            {
                calculer_u_exact(u_e, u_exact, k);
                erreur_infty_k = norme_infty_diff(u_exact, u, nb_pt * nb_pt);
                if (erreur_infty_k > erreur_infty){
                    erreur_infty = erreur_infty_k;
                }
            }
            permut = u; u = u_anc; u_anc = permut;

        }

        terminaison(&permut, &u, &u_anc);

    }

    free(u_exact);

    return erreur_infty;

}