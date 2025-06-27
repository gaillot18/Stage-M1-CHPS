# include <stdio.h>
# include <stdlib.h>
# include <sys/time.h>
# include <mpi.h>
# include <math.h>

# include "../../Librairies/parallele-2.h"

// ======================================================
// Déclarations des variables globales
// ======================================================
// MPI
int rang;
int nb_cpu;
int nb_pt_div_i;
int nb_pt_div_j;
int i_debut;
int i_fin;
int j_debut;
int j_fin;
// Communicateur 2D
MPI_Comm comm_2D;
int dims[2];
int coords[2];
int voisins[4];
int nb_bord_libre;
MPI_Datatype ligne;
MPI_Datatype colonne;
MPI_Datatype bloc_send;
int etiquette = 1;
MPI_Status statut;
// Variable égale pour chaque rang
int N;
int nb_pt;
int nb_iteration;



int main(int argc, char **argv){

    // ======================================================
    // Déclarations des variables
    // ======================================================
    // Temps
    double temps_debut;
    double temps_fin;
    double temps;
    // Buffers MPI
    double *f_div;
    double *u_div;
    // Buffers rang 0
    double *u = NULL;
    double *u_exact = NULL;
    // Autres résultats
    double erreur_infty;
    // Fichiers
    const char *entete;
    double resultats[6];
    const char *nom_fichier_txt;
    // Paramètres
    if (argc > 1){
        N = atoi(argv[1]);
    }
    else{
        N = 10;
    }
    nb_pt = N + 1;


    // ======================================================
    // Initialisation de MPI
    // ======================================================
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rang);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_cpu);
    if (rang == 0){
        printf("------------------------------------------------------------\n");
        printf("Exécution parallèle (pour %d processus) de : parallele-2 (version 3 - méthode itérative - MPI bloquant)\n", nb_cpu);
    }



    // ======================================================
    // Récupération des informations de chaque processus, création de la topologie et des types
    // ======================================================
    creer_topologie();
    infos_topologie();
    infos_processus();
    creer_types();
    if (rang == 0){
        printf("Informations de chaque processus :\n");
    }
    printf("rang = %d, i_debut = %d, j_debut = %d, nb_pt_div_i = %d, nb_pt_div_j = %d, voisins[gauche] = %d, voisins[haut] = %d, "
        "voisins[droite] = %d, voisins[bas] = %d\n", rang, i_debut, j_debut, nb_pt_div_i, nb_pt_div_j, voisins[0], voisins[1],
        voisins[2], voisins[3]);



    // ======================================================
    // Calcul de f et u_exact
    // ======================================================
    f_1(&f_div);
    u_div = (double *)malloc((nb_pt_div_i + 2) * (nb_pt_div_j + 2) * sizeof(double));
    if (rang == 0){
        u = (double *)malloc(nb_pt * nb_pt * sizeof(double));
        u_exact = (double *)malloc(nb_pt * nb_pt * sizeof(double));
        calculer_u_exact(u_e_1, u_exact);
    }
    


    // ======================================================
    // Calcul de u avec mesure de temps
    // ======================================================
    temps_debut = MPI_Wtime();
    calculer_u_jacobi(f_div, u_div);
    regrouper_u(u_div, u);
    temps_fin = MPI_Wtime();
    temps = temps_fin - temps_debut;



    // ======================================================
    // Affichage des informations du problème et des résultats
    // ======================================================
    if (rang == 0){
        erreur_infty = norme_infty_diff(u, u_exact, nb_pt * nb_pt); printf("\n");
        printf("N = %d\nnb_pt * nb_pt = %d\nnb_iteration = %d, erreur_infty = %f\ntemps = %f sec\n", N, nb_pt * nb_pt, nb_iteration, erreur_infty, temps);
    }



    // ======================================================
    // Sauvegarde des résultats dans un fichier
    // ======================================================
    if (rang == 0){
        nom_fichier_txt = "./Textes/resultats.txt";
        entete = "version nb_cpu N nb_iteration erreur_infty temps";
        resultats[0] = 3.0; resultats[1] = (double)nb_cpu; resultats[2] = (double)N; resultats[3] = (double)nb_iteration; resultats[4] = erreur_infty; resultats[5] = temps;
        ecrire_resultats(resultats, entete, 6, nom_fichier_txt);
    }



    // ======================================================
    // Libérations de la mémoire
    // ======================================================
    MPI_Type_free(&ligne);
    MPI_Type_free(&colonne);
    MPI_Type_free(&bloc_send);
    MPI_Comm_free(&comm_2D);
    free(u_div);
    free(f_div);
    if (rang == 0){
        free(u_exact);
        free(u);
    }



    // ======================================================
    // Fermeture de MPI
    // ======================================================
    if (rang == 0){
        printf("Exécution terminée\n");
        printf("------------------------------------------------------------\n");
    }
    MPI_Finalize();
    
    return 0;
    
}