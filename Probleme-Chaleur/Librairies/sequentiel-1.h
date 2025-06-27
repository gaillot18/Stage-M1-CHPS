// ======================================================
// Variables globales
// ======================================================
extern FILE *descripteur;
extern const char *nom_fichier_bin;
extern double L;
extern int N;
extern double h;
extern double T;
extern int N_t;
extern double h_t;
extern int nb_pt;
extern double a;
extern double alpha;
extern double beta;
extern double lambda;



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
void ecrire_double(char *nom_fichier_data, double *t, int n);
void ecrire_resultats(double *resultats, const char *entete, int n, const char *nom_fichier);



// ======================================================
// ../Source/sequentiel-1/resolution.c
// ======================================================
double u_1(double x, double y, double t);
void calculer_u(double *u);
double calculer_u_u_exact(double *u);