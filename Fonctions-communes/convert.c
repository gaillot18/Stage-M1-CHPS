# include <stdio.h>
# include <stdlib.h>
# ifdef USE_MPI
# include <mpi.h>
# endif



// Convertir un fichier data contenant un tableau de doubles en fichier txt
void convertir_data_vers_txt(const char *nom_fichier_data, const char *nom_fichier_txt){

    FILE *fichier_data = fopen(nom_fichier_data, "rb");
    FILE *fichier_txt = fopen(nom_fichier_txt, "w");

    if (fichier_data == NULL || fichier_txt == NULL){
        perror("Erreur d'ouverture du fichier");
    }

    double valeur;
    while (fread(&valeur, sizeof(double), 1, fichier_data) == 1){
        fprintf(fichier_txt, "%f\n", valeur);
    }

    fclose(fichier_data);
    fclose(fichier_txt);

}



// Écrire un tableau de doubles dans un fichier data
void ecrire_double(char *nom_fichier_data, double *t, int n){

    FILE *descripteur = fopen(nom_fichier_data, "ab");
    fwrite(t, sizeof(double), n, descripteur);
    fclose(descripteur);

}



// Écrire les résultats dans un fichier txt
void ecrire_resultats(double *resultats, const char *entete, int n, const char *nom_fichier){

    FILE *descripteur = fopen(nom_fichier, "a+");
    if (descripteur == NULL) {
        printf("Erreur : impossible d'ouvrir le fichier.\n");
        return;
    }

    fseek(descripteur, 0, SEEK_END);
    long taille = ftell(descripteur);
    if (taille == 0){
        const char *p = entete;
        int colonne = 0;
        while (*p && colonne < n){
            char nom[64] = {0};
            int i = 0;
            while (*p && *p != ' ' && i < 63){
                nom[i ++] = *p ++;
            }
            nom[i] = '\0';

            fprintf(descripteur, "%20s", nom);

            while (*p == ' ') p++;

            colonne ++;
        }
        fprintf(descripteur, "\n");
    }

    for (int i = 0 ; i < n ; i ++){
        fprintf(descripteur, "%20.8f", resultats[i]);
    }
    fprintf(descripteur, "\n");

    fclose(descripteur);

}

