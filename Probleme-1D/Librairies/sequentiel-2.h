// ======================================================
// Variables globales
// ======================================================
extern int N;
extern int nb_pt;



// ======================================================
// Structures
// ======================================================
struct mat_2bandes{
    int N;
    double *diag; // taille N - 1
    double *sous_diag; // taille N - 2
};



// ======================================================
// ../../Fonctions-communes/affichage.c
// ======================================================
void afficher_matrice_carre_double(double *A, int n);
void afficher_matrice_carre_int(int *A, int n);
void afficher_vecteur_double(double *v, int n);
void afficher_vecteur_int(int *v, int n);



// ======================================================
// ../../Fonctions-communes/calcul_mat.c
// ======================================================
void init_matrice_carre_zero(int N, double *A);
void somme_matrice_carre(double alpha, double *A, double beta, double *B, int N, double *C);
void produit_matrice_carre(double alpha, double *A, double *B, int N, double *C);
double norme_L2_diff(double *u, double *v, int n);
double norme_infty_diff(double *u, double *v, int n);
double carre_norme_L2_diff(double *u, double *v, int n);
double norme_L2(double *u, int n);
double carre_norme_L2(double *u, int n);
double norme_infty(double *u, int n);



// ======================================================
// ../../Fonctions-communes/convert.c
// ======================================================
void convertir_data_vers_txt(const char *nom_fichier_data, const char *nom_fichier_txt);
void ecrire_double(char *nom_fichier_data, char *nom_fichier_txt, double *t, int n);
void ecrire_resultats(double *resultats, const char *entete, int n, const char *nom_fichier);



// ======================================================
// ../Source/sequentiel-2/affichage-struct.c
// ======================================================
void afficher_mat_2bandes(struct mat_2bandes *A);
void afficher_mat_2bandes_totale(struct mat_2bandes *A);
void mat_2bandes_vers_mat(struct mat_2bandes *A, double **B);
void mat_2bandes_vers_mat_trans(struct mat_2bandes *A, double **B);



// ======================================================
// ../Source/sequentiel-2/resolution.c
// ======================================================
void f_0(double **f);
double u_e_0(double x);
void f_1(double **f);
double u_e_1(double x);
void init_mat_2bandes(struct mat_2bandes *A);
void liberer_mat_2bandes(struct mat_2bandes *A);
void calculer_u_exact(double (*fonction)(double), double *u);
void calculer_cholesky(struct mat_2bandes *L);
void resoudre_cholesky_descente(struct mat_2bandes *L, double *f, double *y);
void resoudre_cholesky_remontee(struct mat_2bandes *L, double *y, double *u);
void resoudre_cholesky(double *f, double *u);