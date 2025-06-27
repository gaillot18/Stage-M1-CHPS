# include <stdio.h>
# include <stdlib.h>
# include <string.h>
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
int nb_pt_div;
int i_debut;
int i_fin;
MPI_Comm comm_1D;
int dims;
int coords;
int voisins[2];
int nb_bord_libre;
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
    // Informations MPI
    int *nb_elements_recus = NULL;
    int *deplacements = NULL;
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
        printf("Exécution parallèle (pour %d processus) de : parallele-2 (version 3 - méthode itérative - MPI)\n", nb_cpu);
    }

    

    // ======================================================
    // Récupération des informations de chaque processus
    // ======================================================
    creer_topologie();
    infos_topologie();
    infos_processus();
    infos_gather(&deplacements, &nb_elements_recus);
    if (rang == 0){
        printf("Informations de chaque processus :\n");
    }
    printf("rang = %d, i_debut = %d, nb_pt_div = %d, voisins[gauche] = %d, voisins[droite] = %d\n", rang, i_debut, nb_pt_div, voisins[0], voisins[1]);
    


    // ======================================================
    // Calcul de f_divise, u_divise et u_exact
    // ======================================================
    f_1(&f_div);
    u_div = (double *)malloc((nb_pt_div + 2) * sizeof(double));
    if (rang == 0){
        u_exact = (double *)malloc(nb_pt * sizeof(double));
        calculer_u_exact(u_e_1, u_exact);
    }
    


    // ======================================================
    // Calcul de u_divise et u avec mesure de temps
    // ======================================================
    temps_debut = MPI_Wtime();
    calculer_u_jacobi(f_div, u_div);
    if (rang == 0){
        u = (double *)malloc(nb_pt * sizeof(double));
    }
    MPI_Gatherv(&(u_div[1]), nb_pt_div, MPI_DOUBLE, u, nb_elements_recus, deplacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    temps_fin = MPI_Wtime();
    temps = temps_fin - temps_debut;
    


    // ======================================================
    // Affichage des informations du problème et des résultats
    // ======================================================
    if (rang == 0){
        erreur_infty = norme_infty_diff(u, u_exact, nb_pt);
        printf("N = %d\nnb_iteration = %d, erreur_infty = %f\ntemps = %f sec\n", N, nb_iteration, erreur_infty, temps);
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
    free(u_div);
    free(f_div);
    if (rang == 0){
        free(u);
        free(u_exact);
        free(nb_elements_recus);
        free(deplacements);
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