#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

int min(int a, int b)
{
    return a < b ? a : b;
}

int compute_nnz(int n_hole, double h, GridParameters *grid, OutputData *output, PhysicalParameters *physical) 
/* 
 * But
 * ===
 * Calculer le nombre d'éléments non nuls de la matrice A
 *
 * Arguments
 * =========
 * n_hole (input) - nombre d'inconnues dans le trou
 * h (input) - pas de discrétisation selon x
 * grid (input) - paramètres de la grille
 * output (input/output) - données de sortie
 * physical (input) - paramètres physiques
 */
{
    printf("\n------ utils.c -> compute_nnz() ------\n");
    int i_start_x = output->i_start_x;
    int i_end_x = output->i_end_x;
    int i_start_y = output->i_start_y;
    int i_end_y = output->i_end_y;
    int mx = grid->mx;                     /* Nombre de points par direction dans la grille selon x */
    int my = grid->my;                     /* Nombre de points par direction dans la grille selon y */
    int nx = grid->nx;                     /* Nombre de points par direction dans la grille selon x (hors bords) */
    int ny = grid->ny;                     /* Nombre de points par direction dans la grille selon y (hors bords) */

    int n_edge, nnz_hole;
    int nx_hole = i_end_x - i_start_x + 1; /* Nombre de points en x dans le trou */
    int ny_hole = i_end_y - i_start_y + 1; /* Nombre de points en y dans le trou */
    int perimeter = 2 * (nx_hole + ny_hole); /* Périmètre du rectangle entourant le trou (nombre de points sur le bord du trou) */
    printf("- Nombre de points dans le trou -> selon x : nx_hole = %d; selon y : ny_hole = %d\n", nx_hole, ny_hole);

    /* Calcul du nombre d'inconnues qui interagissent avec le bord du domaine */
    n_edge = 0;

    /* Vérifie si les inconnues du bord ouest du trou interagissent avec le bord ouest du domaine */
    if (i_start_x == 0)
        n_edge += ny_hole; /* Ajoute les interactions sur le bord ouest du trou */

    /* Vérifie si les inconnues du bord est du trou interagissent avec le bord est du domaine */
    if (i_end_x  == mx - 2)
        n_edge += ny_hole; /* Ajoute les interactions sur le bord est du trou */

    /* Vérifie si les inconnues du bord sud du trou interagissent avec le bord sud du domaine */
    if (i_start_y == 0)
        n_edge += nx_hole; /* Ajoute les interactions sur le bord sud du trou */

    /* Vérifie si les inconnues du bord nord du trou interagissent avec le bord nord du domaine */
    if (i_end_y == my - 2)
        n_edge += nx_hole; /* Ajoute les interactions sur le bord nord du trou */

    printf("- Nombre de points du trou qui se situe à côté d'un bord du domaine -> n_edge = %d\n", n_edge);
    printf("- Nombre de points sur le perimètre du trou -> perimeter = %d\n", perimeter);

    /* Calcul du nombre d'éléments non nuls (nnz) associés au trou */
    nnz_hole = 5 * n_hole - n_edge;
    /* Chaque point à l'intérieur du trou aurait normalement 5 interactions (coefficients non nuls),
      mais les points sur le bord du trou adjacents au bord du domaine en ont moins. */

    /* Calcul du nombre total d'éléments non nuls dans la matrice globale */
    output->nnz = (5 * nx * ny - 2 * (nx + ny)) - nnz_hole - perimeter;
    /* On soustrait du total :
      - les interactions à l'intérieur du trou (nnz_hole),
      - le périmètre du trou (interactions entre le trou et le reste du domaine). */
    return 0;
}

int compute_grid_values(GridParameters *grid, OutputData *output, PhysicalParameters *physical)
/* 
 * But
 * ===
 * Calculer les valeurs de la grille et les stocker dans un tableau à 2 dimensions
 *
 * Arguments
 * =========
 * grid (input) - paramètres de la grille
 * output (input/output) - données de sortie
 * physical (input) - paramètres physiques
 */
{
    printf("\n------ utils.c -> compute_grid_values() ------\n");
    int mx = grid->mx;
    int my = grid->my;
    FILE *file = fopen("test.txt", "w");

    int i_start_x = output->i_start_x;
    int i_end_x = output->i_end_x;
    int i_start_y = output->i_start_y;
    int i_end_y = output->i_end_y;

    int i, j, ind;

    /* Allocation de la mémoire pour le tableau des valeurs de la grille */
    output->grid_values = (double **)malloc(mx * sizeof(double *));
    for (int i = 0; i < mx; i++)
    {
        output->grid_values[i] = (double *)malloc(my * sizeof(double));
    }
    if (output->grid_values == NULL)
    {
        fprintf(stderr, "Erreur d'allocation de mémoire pour le tableau des valeurs de la grille\n");
        return 1;
    }
    ind = 0;
    for (j=0; j<my; j++)
    {
        for (i=0; i<mx; i++)
        {
            if (i==0 || i==mx-1 || j==0 || j==my-1) /* Bords du domaine */
                output->grid_values[i][j] = physical->ue; /* Tension aux bords du domaine */
            else if (i >= i_start_x + 1 && i <= i_end_x + 1 && j >= i_start_y + 1 && j <= i_end_y + 1) /* Trou */
                output->grid_values[i][j] = physical->ui; /* Tension dans le trou */
            else
            {
                output->grid_values[i][j] = output->x[ind]; /* Valeur de la solution */
                ind++; /* Incrémente l'indice de la solution */
            }
        }
    }

    ind = 0;
    for (j=0; j<my; j++)
    {
        for (i=0; i<mx; i++)
        {
            fprintf(file, "i = %d, j = %d, INDEX = %d, grid = %f\n", i, j, ind, output->grid_values[i][j]);
            ind++;
        }
        fprintf(file, "\n");
    }
    fclose(file);
    return 0;
}


void find_extremes_grid_values(OutputData *output, GridParameters *grid)
/* But
* ===
    * Trouver les valeurs minimales et maximales de la grille
* Arguments
* =========
* output (input) - données calculée par le programme sortie (contient le tableau des valeurs de la grille)
* grid (input) - données relatives à la grille
*/
{
    printf("\n------ utils.c -> find_extremes_grid_values() ------\n");
    int mx = grid->mx;
    int my = grid->my;
    int i, j;
    double max = output->grid_values[0][0];
    double min = output->grid_values[0][0];

    for (i=0; i<mx; i++)
    {
        for (j=0; j<my; j++)
        {
            if (output->grid_values[i][j] < min)
                min = output->grid_values[i][j];
            if (output->grid_values[i][j] > max)
                max = output->grid_values[i][j];
        }
    }
    output->min = min;
    output->max = max;
    printf("- Valeur minimale de la grille : %f\n", min);
    printf("- Valeur maximale de la grille : %f\n", max);
}
