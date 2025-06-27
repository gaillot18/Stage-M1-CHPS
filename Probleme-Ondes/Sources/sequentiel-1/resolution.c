# include <stdio.h>
# include <stdlib.h>
# include <sys/time.h>
# include <math.h>
# include <float.h>
# include <limits.h>

# include "../../Librairies/sequentiel-1.h"

# define pi 3.14159265358979323846
# define IDX(i, j) ((j) * (nb_pt) + (i))

double const_1;



// u_0 dont on connait la solution exacte
double u_0(double x){

    double res = sin(pi * x);

    return res;

}



// u_1 dont on connait la solution exacte
double u_1(double x){

    double res = 0.0;

    return res;

}



// Solution exacte
double u_e(double x, double t){

    double res = cos(pi * t) * sin(pi * x);

    return res;
    
}



// Calculer la solution exacte
void calculer_u_exact(double (*fonction)(double, double), double *u, int k){

    for (int i = 0 ; i < nb_pt ; i ++){
        u[i] = fonction(i * h, k * h_t);
    }

}



// Allouer et initialiser u_0
void init_u_0(double (*fonction)(double), double **u_anc_1){

    *u_anc_1 = (double *)malloc(nb_pt * sizeof(double));

    (*u_anc_1)[0] = 0.0; (*u_anc_1)[nb_pt - 1] = 0.0;

    for (int i = 1 ; i < nb_pt - 1; i ++){
        (*u_anc_1)[i] = fonction(i * h);
    }

}



// Allouer et initialiser u_1 (k = 0)
void init_u_1(double (*fonction)(double), double *u_anc_1, double **u_anc_0){
    
    *u_anc_0 = (double *)malloc(nb_pt * sizeof(double));

    (*u_anc_0)[0] = 0.0; (*u_anc_0)[nb_pt - 1] = 0.0;

    for (int i = 1; i < nb_pt - 1; i++) {
        (*u_anc_0)[i] = h_t * fonction(i * h) + u_anc_1[i];
    }

}




// Appliquer le schéma à un point (k > 0)
static inline __attribute__((always_inline)) double schema(double *u_anc_0, double *u_anc_1, int i, int k){

    // const_1 = pow(c, 2) * pow(h_t, 2) / pow(h, 2)
    double res = -u_anc_1[i] + 2 * (1 - const_1) * u_anc_0[i] + const_1 * (u_anc_0[i - 1] + u_anc_0[i + 1]);

    return res;

}



// Écrire dans un fichier
static inline __attribute__((always_inline, unused)) void ecrire_double_iteration(double *u){

    fwrite(u, sizeof(double), nb_pt, descripteur);

}



// Terminer
void terminaison(double **permut, double **u, double **u_anc_0, double **u_anc_1){

    int nb_permut = 0;

    if (N_t % 3 == 1){
        nb_permut = 2;
    }
    else if (N_t % 3 == 2){
        nb_permut = 1;
    }

    for (int i = 0 ; i < nb_permut ; i ++){
        *permut = *u_anc_1; *u_anc_1 = *u_anc_0; *u_anc_0 = *u; *u = *permut;
    }

    free(*u_anc_0); free(*u_anc_1);

}



// Fonction principale (calcul uniquement de la solution approchée)
void calculer_u(double *u){

    const_1 = pow(c, 2) * pow(h_t, 2) / pow(h, 2);
    double *u_anc_0; double *u_anc_1; double *permut;
    init_u_0(u_0, &u_anc_1); init_u_1(u_1, u_anc_1, &u_anc_0);
    for (int i = 0 ; i < nb_pt ; i ++){
        u[i] = 0.0;
    }

    for (int k = 1 ; k <= N_t ; k ++){

        # ifdef ECRITURE
        ecrire_double_iteration(u_anc_0);
        # endif

        for (int i = 1 ; i < nb_pt - 1 ; i ++){
            u[i] = schema(u_anc_0, u_anc_1, i, k);
        }

        permut = u_anc_1; u_anc_1 = u_anc_0; u_anc_0 = u; u = permut;

    }

    # ifdef ECRITURE
    ecrire_double_iteration(u);
    # endif

    terminaison(&permut, &u, &u_anc_0, &u_anc_1);

}



// Fonction principale (calcul de la solution exacte et de la solution approchée pour avoir l'erreur à chaque itération)
double calculer_u_u_exact(double *u){

    double *u_exact = (double *)malloc(nb_pt * sizeof(double));
    double erreur_infty_k; double erreur_infty = 0.0;

    const_1 = pow(c, 2) * pow(h_t, 2) / pow(h, 2);
    double *u_anc_0; double *u_anc_1; double *permut;
    init_u_0(u_0, &u_anc_1); init_u_1(u_1, u_anc_1, &u_anc_0);
    for (int i = 0 ; i < nb_pt ; i ++){
        u[i] = 0.0;
    }

    for (int k = 1 ; k <= N_t ; k ++){

        for (int i = 1 ; i < nb_pt - 1 ; i ++){
            u[i] = schema(u_anc_0, u_anc_1, i, k);
        }
        
        calculer_u_exact(u_e, u_exact, k);
        erreur_infty_k = norme_infty_diff(u_exact, u, nb_pt);
        if (erreur_infty_k > erreur_infty){
            erreur_infty = erreur_infty_k;
        }

        permut = u_anc_1; u_anc_1 = u_anc_0; u_anc_0 = u; u = permut;

    }

    terminaison(&permut, &u, &u_anc_0, &u_anc_1);
    free(u_exact);
    
    return erreur_infty;

}