# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <sys/time.h>
# include <math.h>
# include <mpi.h>
# include <float.h>
# include <limits.h>

# include "../../Librairies/parallele-2.h"

# define SORTIE

# define pi 3.14159265358979323846

double h_carre;



// f dont on connait la solution exacte (1)
void f_0(double **f){

    *f = (double *)malloc((nb_pt_div + 2) * sizeof(double));
    for (int i = 1 ; i < nb_pt_div + 1 ; i ++){
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

    *f = (double *)malloc((nb_pt_div + 2) * sizeof(double));
    double h = 1.0 / N;
    int i_reel = i_debut;
    
    for (int i = 1 ; i < nb_pt_div + 1 ; i ++){
        (*f)[i] = pi * pi * sin(pi * i_reel * h);
        i_reel ++;
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



// Initialiser u_div_anc
void init_u_div_anc(double **u_div_anc){

    *u_div_anc = (double *)malloc((nb_pt_div + 2) * sizeof(double));
    
    for (int i = 0 ; i < nb_pt_div + 2 ; i ++){
        (*u_div_anc)[i] = 0.0;
    }

}



// Appliquer le schéma à un point
static inline __attribute__((always_inline)) double schema(double *f_div, double *u_div_anc, int i){

    double res = 0.5 * ((u_div_anc[i - 1] + u_div_anc[i + 1]) + h_carre * f_div[i]);

    return res;

}



// Calculer la norme infinie relative
static inline __attribute__((always_inline)) double norme_infty_iteration(double *u_div, double *u_div_anc){

    double norme_nume_div = 0.0;
    double norme_deno_div = 0.0;
    double norme_nume;
    double norme_deno;
    double norme;

    for (int i = 1 ; i < nb_pt_div + 1 ; i ++){
        double diff = fabs(u_div[i] - u_div_anc[i]);
        if (diff > norme_nume_div){
            norme_nume_div = diff;
        }
        if (fabs(u_div_anc[i]) > norme_deno_div){
            norme_deno_div = fabs(u_div_anc[i]);
        }
    }

    MPI_Allreduce(&norme_nume_div, &norme_nume, 1, MPI_DOUBLE, MPI_MAX, comm_1D);
    MPI_Allreduce(&norme_deno_div, &norme_deno, 1, MPI_DOUBLE, MPI_MAX, comm_1D);
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
    int i_boucle_debut, i_boucle_fin;
    double *u_div_anc; double *permut;

    // Vecteur de départ
    init_u_div_anc(&u_div_anc);
    for (int i = 0 ; i < nb_pt_div + 2 ; i ++){
        u_div[i] = 0.0;
    }

    // Bornes des boucles
    infos_bornes_boucles(&i_boucle_debut, &i_boucle_fin);

    // Itérations
    for (int iteration = 0 ; iteration < nb_iteration_max && norme > 1e-10 ; iteration ++){

        // Communication
        echanger_halos(u_div_anc);

        for (int i = i_boucle_debut ; i < i_boucle_fin ; i ++){
            u_div[i] = schema(f_div, u_div_anc, i);
        }

        // Test d'arrêt
        norme = norme_infty_iteration(u_div, u_div_anc);

        permut = u_div; u_div = u_div_anc; u_div_anc = permut; nb_iteration ++;

    }

    terminaison(&permut, &u_div, &u_div_anc);

}