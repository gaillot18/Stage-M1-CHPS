# include <stdio.h>
# include <stdlib.h>
# include <sys/time.h>
# include <math.h>

# include "../../Librairies/base.h"

// ======================================================
// Déclarations des variables globales
// ======================================================
int N;
int nb_pt;



int main(int argc, char **argv){

    // ======================================================
    // Déclarations des variables
    // ======================================================
    // Temps
    struct timeval temps_debut;
    struct timeval temps_fin;
    double temps;
    // Buffers
    double *A;
    double *f;
    double *u;
    double *u_exact;
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
    // Initialisation
    // ======================================================
    printf("------------------------------------------------------------\n");
    printf("Exécution séquentielle de : base (version 0 - méthode directe - séquentiel)\n");



    // ======================================================
    // Calcul de f et u_exact
    // ======================================================
    f_1(&f);
    u = (double *)malloc(nb_pt * sizeof(double));
    u_exact = (double *)malloc(nb_pt * sizeof(double));
    calculer_u_exact(u_e_1, u_exact);



    // ======================================================
    // Calcul de u avec mesure de temps
    // ======================================================
    gettimeofday(&temps_debut, NULL);
    A = construire_matrice();
    resoudre_gauss(A, f, u);
    gettimeofday(&temps_fin, NULL);
    temps = (temps_fin.tv_sec - temps_debut.tv_sec) + (temps_fin.tv_usec - temps_debut.tv_usec) / (double)1000000;



    // ======================================================
    // Affichage des informations du problème et des résultats
    // ======================================================
    erreur_infty = norme_infty_diff(u, u_exact, nb_pt);
    printf("N = %d\nnb_pt * nb_pt = %d\nerreur_infty = %f\ntemps = %f sec\n", N, nb_pt * nb_pt, erreur_infty, temps);



    // ======================================================
    // Sauvegarde des résultats dans un fichier
    // ======================================================
    nom_fichier_txt = "./Textes/resultats.txt";
    entete = "version nb_cpu N nb_iteration erreur_infty temps";
    resultats[0] = 0.0; resultats[1] = NAN; resultats[2] = (double)N; resultats[3] = -1.0; resultats[4] = erreur_infty; resultats[5] = temps;
    ecrire_resultats(resultats, entete, 6, nom_fichier_txt);



    // ======================================================
    // Libérations de la mémoire
    // ======================================================
    free(f);
    free(u_exact);
    free(u);
    free(A);



    // ======================================================
    // Fermeture
    // ======================================================
    printf("Exécution terminée\n");
    printf("------------------------------------------------------------\n");
    
    return 0;
    
}