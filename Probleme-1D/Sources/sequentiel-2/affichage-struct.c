# include <stdio.h>
# include <stdlib.h>

# include "../../Librairies/sequentiel-2.h"

# define IDX(i, j) ((i) * ((N) - 1) + (j))
# define idx_max ((N) - 1)



// Afficher une structure mat_2bandes (compressée)
void afficher_mat_2bandes(struct mat_2bandes *A){
    
    printf("N = %d\n", A -> N);

    int N = A -> N;

    printf("diag      =");
    for (int i = 0 ; i < idx_max ; i ++){
        printf("%10.6f ", (A -> diag)[i]);
    }
    printf("\n");

    printf("sous_diag =");
    for (int i = 0 ; i < idx_max - 1 ; i ++){
        printf("%10.6f ", (A -> sous_diag)[i]);
    }
    printf("\n");

}



// Afficher une strucutre mat_2bandes (décompressée) (esthétique : pas utile pour les calculs)
void afficher_mat_2bandes_totale(struct mat_2bandes *A){

    int N = A -> N;
    double zero = 0.0;

    for (int i = 0 ; i < idx_max ; i ++){
        for (int j = 0 ; j < idx_max ; j ++){
            if (i == j){
                printf("%10.6f ", (A -> diag)[i]);
            }
            else if (i == j + 1){
                printf("%10.6f ", (A -> sous_diag)[j]);
            }
            else{
                printf("%10.6f ", zero);
            }
        }
        printf("\n");
    }
    printf("\n");

}



// Convertir une structure mat_2bandes en matrice carré
void mat_2bandes_vers_mat(struct mat_2bandes *A, double **B){

    int N = A -> N;
    double zero = 0.0;

    *B = (double *)malloc(idx_max * idx_max * sizeof(double));
    for (int i = 0 ; i < idx_max ; i ++){
        for (int j = 0 ; j < idx_max; j ++){
            if (i == j){
                (*B)[IDX(i, j)] = (A -> diag)[i];
            }
            else if (i == j + 1){
                (*B)[IDX(i, j)] = (A -> sous_diag)[j];
            }
            else{
                (*B)[IDX(i, j)] = zero;
            }
        }
    }

}



// Convertir une structure mat_2bandes en matrice carré (transposée)
void mat_2bandes_vers_mat_trans(struct mat_2bandes *A, double **B){

    int N = A -> N;
    double zero = 0.0;

    *B = (double *)malloc(idx_max * idx_max * sizeof(double));
    for (int i = 0 ; i < idx_max ; i ++){
        for (int j = 0 ; j < idx_max ; j ++){
            if (i == j){
                (*B)[IDX(i, j)] = (A -> diag)[i];
            }
            else if (i == j + 1){
                (*B)[IDX(j, i)] = (A -> sous_diag)[j];
            }
            else{
                (*B)[IDX(i, j)] = zero;
            }
        }
    }

}