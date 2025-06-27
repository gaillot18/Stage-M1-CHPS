# include <math.h>
# include <stdlib.h>
# ifdef USE_MPI
# include <mpi.h>
# endif



// Initialiser une matrice carré à 0
void init_matrice_carre_zero(int n, double *A){

    for (int i = 0 ; i < n * n ; i ++){
        A[i] = 0;
    }

}



// Somme de matrices
void somme_matrice_carre(double alpha, double *A, double beta, double *B, int n, double *C){

    for (int i = 0 ; i < n ; i ++){
        for (int j = 0 ; j < n ; j ++){
            C[i * n + j] = alpha * A[i * n + j] + beta * B[i * n + j];
        }
    }

}



// Norme L2 de la différence de deux vecteurs
double norme_L2_diff(double *u, double *v, int n){

    double res = 0;

    for (int i = 0 ; i < n ; i ++){
        double diff = u[i] - v[i];
        res += diff * diff;
    }
    
    res = sqrt(res);

    return res;

}



// Norme infinie de la différence de deux vecteurs
double norme_infty_diff(double *u, double *v, int n){

    double res = 0;

    for (int i = 0 ; i < n ; i ++){
        double diff = fabs(u[i] - v[i]);
        if (diff > res){
            res = diff;
        }
    }

    return res;

}



// Carré de la norme L2 de la différence de deux vecteurs
double carre_norme_L2_diff(double *u, double *v, int n){

    double res = 0;

    for (int i = 0 ; i < n ; i ++){
        double diff = u[i] - v[i];
        res += diff * diff;
    }

    return res;

}



// Norme L2 d'un vecteur
double norme_L2(double *u, int n){

    double res = 0;

    for (int i = 0 ; i < n ; i ++){
        res += u[i] * u[i];
    }
    
    res = sqrt(res);

    return res;

}



// Carré de la norme L2 d'un vecteur
double carre_norme_L2(double *u, int n){

    double res = 0;

    for (int i = 0 ; i < n ; i ++){
        res += u[i] * u[i];
    }

    return res;

}



// Norme infinie d'un vecteur
double norme_infty(double *u, int n){

    double res = 0;

    for (int i = 0 ; i < n ; i ++){
        if (u[i] > res){
            res = u[i];
        }
    }

    return res;

}



// Extraire l'intérieur d'une matrice
void extraire_interieur(double *A, double *A_int, int n){

    int idx = 0;
    for (int i = 1 ; i < n - 1 ; i ++){
        for (int j = 1 ; j < n - 1 ; j ++){
            A_int[idx] = A[i * n + j];
            idx ++;
        }
    }

}



// Mettre l'intérieur d'une matrice
void inserer_interieur(double *A_int, double *A, int n){

    int idx = 0;
    for (int i = 1 ; i < n - 1 ; i ++){
        for (int j = 1 ; j < n - 1 ; j ++){
            A[i * n + j] = A_int[idx];
            idx ++;
        }
    }

}