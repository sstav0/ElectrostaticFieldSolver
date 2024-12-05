#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "compute_density.h"

int initialize_density(double ***rho, GridParameters *grid)
/* 
* But
* ===
    * Initialise le tableau de densité de charge rho à 0
* Arguments
* =========
* rho (output) - tableau de taille nx x ny contenant la densité de charge pour chaque inconnue de la grille
* grid (input) - paramètres de la grille
*/
{
    printf("\n------ compute_density.c -> initialize_density() ------\n");
    int nx = grid->nx;
    int ny = grid->ny;

    /* initialisation de rho */
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
            (*rho)[i][j] = 0.0;
    }
    return 0;
}

int compute_density_disk(double ***rho, double Q, GridParameters *grid, PhysicalParameters *physical)
/* But
 * ===
     * Calculer la densité de charge rho dans un disque de rayon r et de centre (rho_0x, rho_0y)
     * rho = 5 si (x - rho_0x)^2 + (y - rho_0y)^2 <= r^2
* Arguments
* =========
* rho (output) - tableau de taille nx x ny contenant la densité de charge pour chaque inconnue de la grille
* Q (intput) - Charge du disque en Coulomb
* grid (input) - paramètres de la grille
* physical (input) - paramètres physiques
*/
{
    printf("\n------ compute_density.c -> compute_density_disk() ------\n");
    int nx = grid->nx;
    int ny = grid->ny;
    double h = (physical->le) / (grid->mx - 1);
    double l = (physical->Le) / (grid->my - 1);

    double r = (physical->le)/20;
    double rho_0x = (physical->li[0])/2;
    double rho_0y = (physical->Le)/2;
    double fun;
    int count = 0;

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fun = pow((i*l - rho_0x), 2) + pow((j*h - rho_0y), 2);
            if (fun <= pow(r, 2))
            {
                (*rho)[i][j] = Q;
                count += 1;
            }
            else
                (*rho)[i][j] = 0.0;
        }
    }

    printf("- Number of points in the disk = %d", count);
    return 0;
}
