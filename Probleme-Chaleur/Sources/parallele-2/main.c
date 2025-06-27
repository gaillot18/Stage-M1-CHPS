# include <stdio.h>
# include <stdlib.h>
# include <sys/time.h>
# include <math.h>
# include <mpi.h>

# include "../../Librairies/parallele-2.h"

# define pi 3.14159265358979323846
//# define EXACTE
//# define ECRITURE

// ======================================================
// Déclarations des variables globales
// ======================================================
// MPI
MPI_File descripteur;
MPI_Datatype vue;
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
const char *nom_fichier_bin;
double L;
int N;
double h;
double T;
int N_t;
double h_t;
int nb_pt;
double a;
double alpha;
double beta;
double lambda;



int main(int argc, char **argv){

    // ======================================================
    // Déclarations des variables
    // ======================================================
    // Temps
    # ifndef EXACTE
    double temps_debut;
    double temps_fin;
    # endif
    double temps = NAN;
    // Buffers MPI
    double *u_div;
    // Autres résultats
    double erreur_infty;
    // Fichiers
    nom_fichier_bin = "Textes/parallele-2/resultats.bin";
    const char *entete;
    double resultats[8];
    const char *nom_fichier_txt;
    // Paramètres
    if (atoi(argv[1]) == 0){ // Mode saisie du nombre de points des maillages
        L = atof(argv[2]);
        N = atoi(argv[3]);
        T = atof(argv[4]);
        N_t = atoi(argv[5]);
        h = L / N;
        h_t = T / N_t;
    }
    else{ // Mode saisie des pas
        L = atof(argv[2]);
        h = atof(argv[3]);
        T = atof(argv[4]);
        h_t = atof(argv[5]);
        N = L / h;
        N_t = T / h_t;
    }
    a = 1.0;
    nb_pt = N + 1;
    alpha = 1.0 - (4 * a * h_t / pow(h, 2));
    beta = a * h_t / pow(h, 2);
    lambda = 2 * a * pow(pi, 2);


    // ======================================================
    // Initialisation de MPI
    // ======================================================
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rang);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_cpu);
    if (rang == 0){
        printf("------------------------------------------------------------\n");
        printf("Exécution parallèle (pour %d processus) de : parallele-2 (version 3 - schéma explicite - MPI bloquant)\n", nb_cpu);
    }
    # ifdef ARRET
    if (beta > 0.25){
        if (rang == 0){
            printf("beta = %f > 1/4\n", beta);
            printf("Exécution terminée\n");
            printf("------------------------------------------------------------\n");
        }
        MPI_Finalize();
        return 0;
    }
    # endif



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
    // Calcul de u_exact
    // ======================================================
    u_div = (double *)malloc((nb_pt_div_i + 2) * (nb_pt_div_j + 2) * sizeof(double));
    # ifdef EXACTE
    erreur_infty = calculer_u_u_exact(u_div);
    # else
    erreur_infty = NAN;
    # endif
    


    // ======================================================
    // Calcul de u_div avec mesure de temps
    // ======================================================
    # ifndef EXACTE
    {
        temps_debut = MPI_Wtime();
        # ifdef ECRITURE
        MPI_File_open(comm_2D, nom_fichier_bin, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &descripteur);
        # endif
        calculer_u(u_div);
        # ifdef ECRITURE
        MPI_File_close(&descripteur);
        # endif
        temps_fin = MPI_Wtime();
        temps = temps_fin - temps_debut;
    }
    # endif



    // ======================================================
    // Affichage des informations du problème et des résultats
    // ======================================================
    if (rang == 0){
        printf("L = %f, N = %d, nb_pt = %d, nb_pt * nb_pt = %d, h = %f\n", L, N, nb_pt, nb_pt * nb_pt, h);
        printf("T = %f, N_t = %d, h_t = %f\n", T, N_t, h_t);
        printf("alpha = %f, beta = %f\n", alpha, beta);
        printf("erreur_infty = %f\n", erreur_infty);
        printf("temps = %f\n", temps);
        printf("taille fichier (si écriture) = %f Go\n", sizeof(double) * nb_pt * nb_pt * (N_t + 1) * 1e-9);
    }



    // ======================================================
    // Sauvegarde des résultats dans un fichier
    // ======================================================
    if (rang == 0){
        nom_fichier_txt = "./Textes/resultats.txt";
        entete = "version nb_cpu N h N_t h_t erreur_infty temps";
        resultats[0] = 3.0; resultats[1] = nb_cpu; resultats[2] = (double)N; resultats[3] = h; resultats[4] = (double)N_t; resultats[5] = h_t; resultats[6] = erreur_infty, resultats[7] = temps;
        ecrire_resultats(resultats, entete, 8, nom_fichier_txt);
    }



    // ======================================================
    // Libérations de la mémoire
    // ======================================================
    MPI_Type_free(&ligne);
    MPI_Type_free(&colonne);
    MPI_Type_free(&bloc_send);
    MPI_Comm_free(&comm_2D);
    free(u_div);



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