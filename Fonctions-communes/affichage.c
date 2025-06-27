# include <stdio.h>
# ifdef USE_MPI
# include <mpi.h>
# endif



// Afficher une matrice carré de doubles
void afficher_matrice_carre_double(double *A, int n){
    for (int i = 0 ; i < n ; i ++){
        for (int j = 0 ; j < n ; j ++){
            printf("%10.6f ", A[i * n + j]);
        }
        printf("\n");
    }
}



// Afficher une matrice de doubles
void afficher_matrice_double(double *A, int n, int m){
    for (int i = 0 ; i < n ; i ++){
        for (int j = 0 ; j < m ; j ++){
            printf("%10.6f ", A[i * n + j]);
        }
        printf("\n");
    }
}



// Afficher une matrice carré d'int
void afficher_matrice_carre_int(int *A, int n){
    for (int i = 0 ; i < n ; i ++){
        for (int j = 0 ; j < n ; j ++){
            printf("%d ", A[i * n + j]);
        }
        printf("\n");
    }
}



// Afficher une matrice d'int
void afficher_matrice_int(int *A, int n, int m){
    for (int i = 0 ; i < n ; i ++){
        for (int j = 0 ; j < m ; j ++){
            printf("%d ", A[i * n + j]);
        }
        printf("\n");
    }
}



// Afficher un vecteur de doubles
void afficher_vecteur_double(double *v, int n){
    for (int i = 0 ; i < n ; i ++){
        printf("%10.6f ", v[i]);
    }
    printf("\n");
}



// Afficher un vecteur d'int
void afficher_vecteur_int(int *v, int n){
    for (int i = 0 ; i < n ; i ++){
        printf("%d ", v[i]);
    }
    printf("\n");
}