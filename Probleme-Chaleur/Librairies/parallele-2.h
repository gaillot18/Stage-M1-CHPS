// ======================================================
// Variables globales
// ======================================================
extern MPI_File descripteur;
extern MPI_Datatype vue;
extern int rang;
extern int nb_cpu;
extern int nb_pt_div_i;
extern int nb_pt_div_j;
extern int i_debut;
extern int i_fin;
extern int j_debut;
extern int j_fin;
extern MPI_Comm comm_2D;
extern int dims[2];
extern int coords[2];
extern int voisins[4];
extern int nb_bord_libre;
extern MPI_Datatype ligne;
extern MPI_Datatype colonne;
extern MPI_Datatype bloc_send;
extern int etiquette;
extern MPI_Status statut;
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
// ../Source/parallele-2/parallele.c
// ======================================================
void affichage_ordonne(double *u_divise, char *message);
void creer_topologie();
void infos_processus();
void infos_topologie();
void creer_types();
void infos_bord();
void echanger_halos(double *u_div);
void infos_bornes_boucles(int *i_boucle_debut, int *j_boucle_debut, int *i_boucle_fin, int *j_boucle_fin);
void regrouper_u(double *u_div, double *u);



// ======================================================
// ../Source/parallele-2/resolution.c
// ======================================================
double u_1(double x, double y, double t);
void calculer_u(double *u);
double calculer_u_u_exact(double *u);