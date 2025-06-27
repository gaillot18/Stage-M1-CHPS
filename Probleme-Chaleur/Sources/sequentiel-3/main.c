# include <stdio.h>
# include <stdlib.h>
# include <sys/time.h>
# include <cholmod.h>
# include <math.h>

# include "../../Librairies/sequentiel-2.h"

# define pi 3.14159265358979323846
# define idx_max (((N) - 1) * ((N) - 1))
# define nb_elements (((N) - 3) * (5 * (N) + 1) + 12)

// ======================================================
// Déclarations des variables globales
// ======================================================
FILE *descripteur;
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
double gamma;
double lambda;
cholmod_common c;



int main(int argc, char **argv){

    // ======================================================
    // Déclarations des variables
    // ======================================================
    // Temps
    # ifndef EXACTE
    struct timeval temps_debut;
    struct timeval temps_fin;
    # endif
    double temps = NAN;
    // Buffers
    double *u;
    // Autres résultats
    double erreur_infty;
    // Fichiers
    nom_fichier_bin = "Textes/sequentiel-3/resultats.bin";
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
    gamma = a * h_t / (2.0 * pow(h, 2));
    alpha = 1.0 + 4.0 * gamma;
    beta = 1.0 - 4.0 * gamma;
    lambda = 2 * a * pow(pi, 2);
    // Cholmod
    int *lignes;
    double *valeurs;
    int *offsets;
    cholmod_sparse *A;


    // ======================================================
    // Initialisation de Cholmod
    // ======================================================
    printf("------------------------------------------------------------\n");
    printf("Exécution séquentielle de : sequentiel-3 (version 5 - schéma semi-implicite - séquentiel)\n");
    cholmod_start(&c);



    // ======================================================
    // Calcul de u_exact
    // ======================================================
    u = (double *)malloc(nb_pt * nb_pt * sizeof(double));
    # ifdef EXACTE
    construire_matrice_creuse(&lignes, &valeurs, &offsets);
    A = init_matrice_creuse(offsets, lignes, valeurs);
    erreur_infty = resoudre_calculer_u_exact(A, u);
    # else
    erreur_infty = NAN;
    # endif



    // ======================================================
    // Calcul de u avec mesure de temps
    // ======================================================
    # ifndef EXACTE
    {
        gettimeofday(&temps_debut, NULL);
        # ifdef ECRITURE
        descripteur = fopen(nom_fichier_bin, "ab");
        # endif
        construire_matrice_creuse(&lignes, &valeurs, &offsets);
        A = init_matrice_creuse(offsets, lignes, valeurs);
        resoudre(A, u);
        # ifdef ECRITURE
        fclose(descripteur);
        # endif
        gettimeofday(&temps_fin, NULL);
        temps = (temps_fin.tv_sec - temps_debut.tv_sec) + (temps_fin.tv_usec - temps_debut.tv_usec) / (double)1000000;
    }
    # endif



    // ======================================================
    // Affichage des informations du problème et des résultats
    // ======================================================
    printf("L = %f, N = %d, nb_pt = %d, nb_pt * nb_pt = %d, h = %f\n", L, N, nb_pt, nb_pt * nb_pt, h);
    printf("T = %f, N_t = %d, h_t = %f\n", T, N_t, h_t);
    printf("alpha = %f, beta = %f, gamma = %f\n", alpha, beta, gamma);
    printf("erreur_infty = %f\n", erreur_infty);
    printf("temps = %f\n", temps);
    printf("taille fichier (si écriture) = %f Go\n", sizeof(double) * nb_pt * nb_pt * (N_t + 1) * 1e-9);

    

    // ======================================================
    // Sauvegarde des résultats dans un fichier
    // ======================================================
    nom_fichier_txt = "./Textes/resultats.txt";
    entete = "version nb_cpu N h N_t h_t erreur_infty temps";
    resultats[0] = 5.0; resultats[1] = NAN; resultats[2] = (double)N; resultats[3] = h; resultats[4] = (double)N_t; resultats[5] = h_t; resultats[6] = erreur_infty, resultats[7] = temps;
    ecrire_resultats(resultats, entete, 8, nom_fichier_txt);



    // ======================================================
    // Libérations de la mémoire
    // ======================================================
    cholmod_free_sparse(&A, &c);
    free(u);



    // ======================================================
    // Fermeture de Cholmod
    // ======================================================
    cholmod_finish(&c);
    printf("Exécution terminée\n");
    printf("------------------------------------------------------------\n");
    
    return 0;
    
}