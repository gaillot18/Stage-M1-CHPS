# include <stdio.h>
# include <stdlib.h>

# include "../../Librairies/sequentiel-2.h"

# define IDX(i, j) ((i) * (((N) - 1) * ((N) - 1)) + (j))
# define idx_max (((N) - 1) * ((N) - 1))



// Afficher une structure mat_Nbandes (compressée)
void afficher_mat_Nbandes(struct mat_Nbandes *A){
    
    printf("N = %d\n", A -> N);

    int N = A -> N;

    for (int d = 0 ; d < N ; d ++){
        printf("diag[%d] = ", d);
        for (int j = 0 ; j < (N - 1) * (N - 1) && j + d < (N - 1) * (N - 1) ; j ++){
            printf("%10.6f ", (A -> diags)[d][j]);
        }
        printf("\n");
    }

}



// Afficher une strucutre mat_2bandes (décompressée) (esthétique : pas utile pour les calculs)
void afficher_mat_Nbandes_totale(struct mat_Nbandes *A){

    int N = A -> N;
    double zero = 0.0;

    for (int i = 0 ; i < idx_max ; i ++){
        for (int j = 0 ; j < idx_max ; j ++) {
            int d = i - j;
            if (d >= 0 && d < N && j < idx_max - d){
                printf("%10.6f ", A -> diags[d][j]);
            }
            else
            {
                printf("%10.6f ", zero);
            }
        }
        printf("\n");
    }
    printf("\n");

}



// Convertir une structure mat_2bandes en matrice carré
void mat_Nbandes_vers_mat(struct mat_Nbandes *A, double **B){

    int N = A -> N;
    *B = (double *)malloc(idx_max * idx_max * sizeof(double));
    double zero = 0.0;

    for (int i = 0 ; i < idx_max ; i ++){
        for (int j = 0 ; j < idx_max ; j ++){
            int d = i - j;
            if (d >= 0 && d < N && j < idx_max - d){
                (*B)[IDX(i, j)] = (A -> diags)[d][j];
            }
            else
            {
                (*B)[IDX(i, j)] = zero;
            }
        }
    }

}



// Convertir une structure mat_2bandes en matrice carré (transposée)
void mat_Nbandes_vers_mat_trans(struct mat_Nbandes *A, double **B){

    int N = A -> N;
    *B = (double *)malloc(idx_max * idx_max * sizeof(double));
    double zero = 0.0;

    for (int i = 0 ; i < idx_max ; i ++){
        for (int j = 0 ; j < idx_max ; j ++) {
            int d = i - j;
            if (d >= 0 && d < N && j < idx_max - d){
                (*B)[IDX(j, i)] = (A -> diags)[d][j];
            }
            else
            {
                (*B)[IDX(i, j)] = zero;
            }
        }
    }

}