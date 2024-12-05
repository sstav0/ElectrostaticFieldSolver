#ifndef UTILS_H
#define UTILS_H

int min(int a, int b);

/* Strucutre reprenant les variables relatives à la grille */
typedef struct {
    int mx; /* Nombre de points en x */
    int my; /* Nombre de points en y */
    int nx; /* Nombre de cellules en x */
    int ny; /* Nombre de cellules en y */
} GridParameters;

/* Structure reprenant les paramètres physiques du problème */
typedef struct {
    int ue; /* Valeur de la tension sur les bords du domaine */
    int ui; /* Valeur de la tension sur les bords du trou et à l'intérieur de celui-ci */
    double le; /* Longueur du domaine */
    double Le; /* Largeur du domaine */
    double li[2]; /* Coordonnées du trou selon x */
    double Li[2]; /* Coordonnées du trou selon y */
    double fac_eps; /* Facteur de permittivité */
    int refinement; /* Raffinement de la grille */
} PhysicalParameters;

/* Structure reprenant les données de sortie du programme */
typedef struct {
    int n; /* Nombre de points de la grille */
    int nnz; /* Nombre de valeurs non nulles dans la matrice A */
    int i_start_x; /* Indice du point qui se situe en bas à gauche du trou */
    int i_end_x;  /* Indice du point qui se situe en bas à droite du trou */
    int i_start_y; /* Indice du point qui se situe en bas à gauche du trou */
    int i_end_y; /* Indice du point qui se situe en haut à droite du trou */
    int *ia; /* Tableau des indices de début de ligne */
    int *ja; /* Tableau des indices de colonne */
    double *a; /* Tableau des valeurs non nulles de la matrice A */
    double *b; /* Tableau des valeurs de la matrice B */
    double *x; /* Tableau des valeurs de la matrice de solutions */
    double **grid_values; /* Tableau des valeurs pour chaque point de la grille */
    double ***EField; /* Tableau des valeurs du champ électrique pour chaque point de la grille */
    double flux; /* Flux du champ électrique à travers un contour englobant le trou */
    double min; /* Valeur minimale du potentiel électrique de la grille */
    double max; /* Valeur maximale du potentiel électrique de la grille */
} OutputData;

int compute_grid_values(GridParameters *grid, OutputData *output, PhysicalParameters *physical);
int compute_nnz(int n_hole, double h, GridParameters *grid, OutputData *output, PhysicalParameters *physical);
void find_extremes_grid_values(OutputData *output, GridParameters *grid);
#endif
