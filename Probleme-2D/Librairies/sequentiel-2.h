// ======================================================
// Variables globales
// ======================================================
extern int N;
extern int nb_pt;



// ======================================================
// Structures
// ======================================================
struct mat_Nbandes{
    int N;
    double **diags;
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
void extraire_interieur(double *A, double *A_int, int n);
void inserer_interieur(double *A_int, double *A, int n);



// ======================================================
// ../../Fonctions-communes/convert.c
// ======================================================
void convertir_data_vers_txt(const char *nom_fichier_data, const char *nom_fichier_txt);
void ecrire_double(char *nom_fichier_data, char *nom_fichier_txt, double *t, int n);
void ecrire_resultats(double *resultats, const char *entete, int n, const char *nom_fichier);



// ======================================================
// ../Source/sequentiel-2/affichage-struct.c
// ======================================================
void afficher_mat_Nbandes(struct mat_Nbandes *A);
void afficher_mat_Nbandes_totale(struct mat_Nbandes *A);
void mat_Nbandes_vers_mat(struct mat_Nbandes *A, double **B);
void mat_Nbandes_vers_mat_trans(struct mat_Nbandes *A, double **B);



// ======================================================
// ../Source/sequentiel-2/resolution.c
// ======================================================
void f_1(double **f);
double u_e_1(double x, double y);
void init_mat_Nbandes(struct mat_Nbandes *A);
void liberer_mat_Nbandes(struct mat_Nbandes *A);
void calculer_cholesky(struct mat_Nbandes *L);
void calculer_u_exact(double (*fonction)(double, double), double *u);
void resoudre_cholesky(double *f, double *u);