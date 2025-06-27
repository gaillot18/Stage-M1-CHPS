# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <sys/time.h>
# include <mpi.h>

# include "../../Librairies/parallele-2.h"



void affichage_ordonne(double *u_div, char *message){
    for (int i = 0 ; i < nb_cpu ; i ++){
            if (rang == i){
                printf("rang = %d, %s :\n", rang, message);
                afficher_vecteur_double(u_div, nb_pt_div);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
}



// Créer la topologie cartésienne 1D du domaine [0, 1]
void creer_topologie(){

    int tore = 0;
    dims = 0;

    MPI_Dims_create(nb_cpu, 1, &dims);

    MPI_Cart_create(MPI_COMM_WORLD, 1, &dims, &tore, 0, &comm_1D);

    MPI_Barrier(comm_1D);

}



// Obtenir des informations sur la topologie
void infos_topologie(){

    MPI_Cart_coords(comm_1D, rang, 1, &coords);

    MPI_Cart_shift(comm_1D, 0, 1, &(voisins[0]), &(voisins[1]));

    nb_bord_libre = 2;
    for (int i = 0 ; i < 2 ; i ++){
        if (voisins[i] == -1){
            nb_bord_libre --;
        }
    }

    MPI_Barrier(comm_1D);

}



// Obtenir les intervalles à traiter de chaque processus
void infos_processus(){

    i_debut = (coords * nb_pt) / dims;
    i_fin = ((coords + 1) * nb_pt) / dims - 1;
    nb_pt_div = i_fin - i_debut + 1;

    MPI_Barrier(comm_1D);

}



void infos_bornes_boucles(int *i_boucle_debut, int *i_boucle_fin){

    *i_boucle_debut = 1;
    *i_boucle_fin = nb_pt_div + 1;

    if (i_debut == 0){
        (*i_boucle_debut) ++;
    }
    if (i_fin == nb_pt - 1){
        (*i_boucle_fin) --;
    }

}



void infos_gather(int **deplacements, int **nb_elements_recus){

    if (rang == 0){

        *deplacements = (int *)malloc(nb_cpu * sizeof(int));
        *nb_elements_recus = (int *)malloc(nb_cpu * sizeof(int));

        for (int i = 0; i < nb_cpu; i++){
            int i_debut_i = (i * nb_pt) / nb_cpu;
            int i_fin_i = ((i + 1) * nb_pt) / nb_cpu - 1;
            int nb_pt_div_i = i_fin_i - i_debut_i + 1;

            (*deplacements)[i] = i_debut_i;
            (*nb_elements_recus)[i] = nb_pt_div_i;
        }

    }

}



// Effectuer les communications des cellules fantômes
void echanger_halos(double *u_div){

    // Envoi gauche, reception droite
    MPI_Sendrecv(&(u_div[1]), 1, MPI_DOUBLE, voisins[0], etiquette, &(u_div[nb_pt_div + 1]), 1, MPI_DOUBLE, voisins[1], etiquette, comm_1D, &statut);

    // Envoi droite, reception gauche
    MPI_Sendrecv(&(u_div[nb_pt_div]), 1, MPI_DOUBLE, voisins[1], etiquette, &(u_div[0]), 1, MPI_DOUBLE, voisins[0], etiquette, comm_1D, &statut);
    
}