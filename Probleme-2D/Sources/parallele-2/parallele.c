# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <sys/time.h>
# include <mpi.h>

# include "../../Librairies/parallele-2.h"

# define IDX(i, j) ((j) * (nb_pt_div_i + 2) + (i))



// (Debug) Afficher un vecteur dans l'ordre des rangs avec un message
void affichage_ordonne(double *u_divise, char *message){
    for (int i = 0 ; i < nb_cpu ; i ++){
        if (rang == i){
            printf("rang = %d, %s :\n", rang, message);
            //afficher_vecteur_double(u_divise, nb_pt_divise);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}



// Créer la topologie cartésienne 2D du domaine [0, 1] x [0, 1]
void creer_topologie(){

    int tore[2] = {0, 0};
    dims[0] = 0; dims[1] = 0;

    MPI_Dims_create(nb_cpu, 2, dims);

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, tore, 0, &comm_2D);

    MPI_Barrier(comm_2D);

}



// Obtenir des informations sur la topologie
void infos_topologie(){

    MPI_Cart_coords(comm_2D, rang, 2, coords);

    MPI_Cart_shift(comm_2D, 0, 1, &(voisins[0]), &(voisins[2]));
    MPI_Cart_shift(comm_2D, 1, 1, &(voisins[3]), &(voisins[1]));

    nb_bord_libre = 4;
    for (int i = 0 ; i < 4 ; i ++){
        if (voisins[i] == -1){
            nb_bord_libre --;
        }
    }

    MPI_Barrier(comm_2D);

}



// Obtenir les intervalles à traiter de chaque processus
void infos_processus(){

    i_debut = (coords[0] * nb_pt) / dims[0];
    i_fin = ((coords[0] + 1) * nb_pt) / dims[0] - 1;
    nb_pt_div_i = i_fin - i_debut + 1;

    j_debut = (coords[1] * nb_pt) / dims[1];
    j_fin = ((coords[1] + 1) * nb_pt) / dims[1] - 1;
    nb_pt_div_j = j_fin - j_debut + 1;

    MPI_Barrier(comm_2D);

}



// Créer les types dérivés lignes, colonnes et bloc_send
void creer_types(){

    int taille_send[2] = {nb_pt_div_j + 2, nb_pt_div_i + 2};
    int sous_taille_send[2] = {nb_pt_div_j, nb_pt_div_i};
    int debut_send[2] = {1, 1};

    MPI_Type_contiguous(nb_pt_div_i, MPI_DOUBLE, &ligne);
    MPI_Type_vector(nb_pt_div_j, 1, nb_pt_div_i + 2, MPI_DOUBLE, &colonne);

    MPI_Type_create_subarray(2, taille_send, sous_taille_send, debut_send, MPI_ORDER_C, MPI_DOUBLE, &bloc_send);

    MPI_Type_commit(&ligne);
    MPI_Type_commit(&colonne);
    MPI_Type_commit(&bloc_send);

    MPI_Barrier(comm_2D);

}



// Obtenir les indices de départ et d'arrivé de la boucle principale du schéma (adaptés pour itérer sur les bords locaux qui ne sont pas globaux)
void infos_bornes_boucles(int *i_boucle_debut, int *j_boucle_debut, int *i_boucle_fin, int *j_boucle_fin){

    *i_boucle_debut = 1; *j_boucle_debut = 1;
    *i_boucle_fin = nb_pt_div_i + 1; *j_boucle_fin = nb_pt_div_j + 1;

    if (i_debut == 0){
        (*i_boucle_debut) ++;
    }
    if (j_debut == 0){
        (*j_boucle_debut) ++;
    }
    if (i_fin == nb_pt - 1){
        (*i_boucle_fin) --;
    }
    if (j_fin == nb_pt - 1){
        (*j_boucle_fin) --;
    }

}



// Effectuer les communications des cellules fantômes
void echanger_halos(double *u_div){

    // Envoi haut, reception bas
    MPI_Sendrecv(&(u_div[IDX(1, nb_pt_div_j)]), 1, ligne, voisins[1], etiquette, &(u_div[IDX(1, 0)]), 1, ligne, voisins[3], etiquette, comm_2D, &statut);

    // Envoi bas, reception haut
    MPI_Sendrecv(&(u_div[IDX(1, 1)]), 1, ligne, voisins[3], etiquette, &(u_div[IDX(1, nb_pt_div_j + 1)]), 1, ligne, voisins[1], etiquette, comm_2D, &statut);

    // Envoi gauche, reception droite
    MPI_Sendrecv(&(u_div[IDX(1, 1)]), 1, colonne, voisins[0], etiquette, &(u_div[IDX(nb_pt_div_i + 1, 1)]), 1, colonne, voisins[2], etiquette, comm_2D, &statut);

    // Envoi droite, reception gauche
    MPI_Sendrecv(&(u_div[IDX(nb_pt_div_i, 1)]), 1, colonne, voisins[2], etiquette, &(u_div[IDX(0, 1)]), 1, colonne, voisins[0], etiquette, comm_2D, &statut);

}



// Regrouper les parties finales dans un vecteur sur le rang 0
// Principe : Le rang 0 crée dynamiquement un type dérivé pour recevoir le domaine adapté à la taille de chaque processus
void regrouper_u(double *u_div, double *u){

    // Dans l'ordre des rangs
    for (int i = 1 ; i < nb_cpu ; i ++){

        if (rang == i){

            // Envoyer les informations nécessaire au rang 0 pour construire le type dérivé dynamiquement
            MPI_Send(&nb_pt_div_i, 1, MPI_INT, 0, etiquette, comm_2D);
            MPI_Send(&nb_pt_div_j, 1, MPI_INT, 0, etiquette, comm_2D);
            MPI_Send(&i_debut, 1, MPI_INT, 0, etiquette, comm_2D);
            MPI_Send(&j_debut, 1, MPI_INT, 0, etiquette, comm_2D);
            MPI_Send(u_div, 1, bloc_send, 0, etiquette, comm_2D);

        }
        else if (rang == 0){

            MPI_Datatype bloc_recv;

            int nb_pt_div_i_recv;
            int nb_pt_div_j_recv;
            int i_debut_recv;
            int j_debut_recv;

            // Recevoir les informations nécessaire des autres rangs pour construire le type dérivé dynamiquement
            MPI_Recv(&nb_pt_div_i_recv, 1, MPI_INT, i, etiquette, comm_2D, &statut);
            MPI_Recv(&nb_pt_div_j_recv, 1, MPI_INT, i, etiquette, comm_2D, &statut);
            MPI_Recv(&i_debut_recv, 1, MPI_INT, i, etiquette, comm_2D, &statut);
            MPI_Recv(&j_debut_recv, 1, MPI_INT, i, etiquette, comm_2D, &statut);

            int taille_recv[2] = {nb_pt, nb_pt};
            int sous_taille_recv[2] = {nb_pt_div_j_recv, nb_pt_div_i_recv};
            int debut_recv[2] = {j_debut_recv, i_debut_recv};

            // Type dérvié dynamique
            MPI_Type_create_subarray(2, taille_recv, sous_taille_recv, debut_recv, MPI_ORDER_C, MPI_DOUBLE, &bloc_recv);
            MPI_Type_commit(&bloc_recv);

            MPI_Recv(u, 1, bloc_recv, i, etiquette, comm_2D, &statut);

            MPI_Type_free(&bloc_recv);

        }

        MPI_Barrier(comm_2D);

    }

    if (rang == 0){

        // Le rang 0 connait déjà ses propres informations, pas d'envois / récéptions

        MPI_Datatype bloc_recv;

        int taille_recv[2] = {nb_pt, nb_pt};
        int sous_taille_recv[2] = {nb_pt_div_j, nb_pt_div_i};
        int debut_recv[2] = {j_debut, i_debut};
        
        MPI_Type_create_subarray(2, taille_recv, sous_taille_recv, debut_recv, MPI_ORDER_C, MPI_DOUBLE, &bloc_recv);
        MPI_Type_commit(&bloc_recv);

        MPI_Sendrecv(u_div, 1, bloc_send, 0, etiquette, u, 1, bloc_recv, 0, etiquette, comm_2D, &statut);

        MPI_Type_free(&bloc_recv);

    }

    MPI_Barrier(comm_2D);

}