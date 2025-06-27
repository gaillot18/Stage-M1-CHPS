// ======================================================
// Variables globales
// ======================================================
extern int rang;
extern int nb_cpu;
extern int nb_pt_div;
extern int i_debut;
extern int i_fin;
extern MPI_Comm comm_1D;
extern int dims;
extern int coords;
extern int voisins[2];
extern int nb_bord_libre;
extern int coins[4];
extern int etiquette;
extern MPI_Status statut;
extern int N;
extern int nb_pt;
extern int nb_iteration;



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
// ../Source/parallele-2/parallele.c
// ======================================================
void affichage_ordonne(double *u_div, char *message);
void creer_topologie();
void infos_topologie();
void infos_processus();
void echanger_halos(double *u_div);
void infos_bornes_boucles(int *i_boucle_debut, int *i_boucle_fin);
void infos_gather(int **deplacements, int **nb_elements_recus);



// ======================================================
// ../Source/parallele-2/resolution.c
// ======================================================
void f_0(double **f);
double u_e_0(double x);
void f_1(double **f);
double u_e_1(double x);
void calculer_u_exact(double (*fonction)(double), double *u);
void calculer_u_jacobi(double *f, double *u);
