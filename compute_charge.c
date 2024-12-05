#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "compute_charge.h"

double compute_EField_x(double u1, double u2, double h)
/* But
 * ===
     * Calculer le champ électrique selon x pour chaque couple de points
 * Arguments
 * =========
 * u1 (input) - valeur du potentiel électrique au point i
 * u2 (input) - valeur du potentiel électrique au point i+1
 * h (input) - pas de discrétisation selon x
 */
{
    double EField_x = -((u2 - u1) / (h*1e-3)); //E_x = - dV/dx
    /*printf("EField_x = -(%f-%f)/%f = %f\n", u2, u1, h*1e-3, EField_x);*/
    return EField_x;
}

double compute_EField_y(double u1, double u2, double l)
/* But
 * ===
     * Calculer le champ électrique selon y pour chaque couple de points
 * Arguments
 * =========
 * u1 (input) - valeur du potentiel électrique au point i
 * u2 (input) - valeur du potentiel électrique au point i+1
 * l (input) - pas de discrétisation selon y
 */
{
    double EField_y = - ((u2 - u1) / (l*1e-3)); //E_x = - dV/dx
    /*printf("EField_y = -(%f-%f)/%f = %f\n", u2, u1, l*1e-3, EField_y);*/
    return EField_y;
}

int compute_EField(GridParameters *grid, PhysicalParameters *physical, OutputData *output)
/* But
 * ===
     * Calculer le champ électrique pour chaque point de la grille
 * Arguments
 * =========
 * grid (input) - paramètres de la grille
 * physical (input) - paramètres physiques
 * output (input) - données de sortie
 */
{
    printf("\n------ compute_charge.c -> compute_EField() ------\n");
    int i, j;
    double h = physical->Le / (grid->my - 1);
    double l = physical->le / (grid->mx - 1);

    /* Allouer de la mémoire pour un tableau à 2D contenant pour chaque point, les coordonées x et y d'un vecteur du champ élec. */
    output->EField = (double ***)malloc((grid->nx + 1) * sizeof(double **)); 
    if (output->EField == NULL)
    {
        printf("Erreur d'allocation de mémoire pour le tableau EField\n");
        return 1;
    }

    for (int i=0; i<=grid->nx; i++) 
    {
        output->EField[i] = (double **)malloc((grid->ny + 1) * sizeof(double*));
        if (output->EField[i] == NULL)
        {
            printf("Erreur d'allocation de mémoire pour le tableau EField\n");
            return 1;
        }
        for (int j=0; j<=grid->ny; j++)
        {
            output->EField[i][j] = (double*)malloc(2*sizeof(double));
            if (output->EField[i][j] == NULL)
            {
                printf("Erreur d'allocation de mémoire pour le tableau EField\n");
                return 1;
            }
        }
    }

    /* Remplissage des valeurs de potentiel électrique */
    if (compute_grid_values(grid, output, physical))
    {
        printf("Erreur lors du calcul des valeurs de la grille\n");
        return 1;
    }

    /* Calcul du champ électrique et stockage des valeurs dans le tableau EField */
    for (i = 0; i <= grid->nx; i++) /* On s'arrête à l'avant dernier point de la grille selon les axes x et y */
    {
        for (j = 0; j <= grid->ny; j++)
        {
            output->EField[i][j][0] = compute_EField_x((output->grid_values)[i][j], (output->grid_values)[i+1][j], h);
            output->EField[i][j][1] = compute_EField_y((output->grid_values)[i][j], (output->grid_values)[i][j+1], l);
        }
    }

    return 0;
}

int compute_flux(int cx[2], int cy[2], GridParameters *grid, PhysicalParameters *physical, OutputData *output)
/* But
 * ===
     * Calculer le flux du champ électrique à travers un contour donné
     * J'utilise l'approximation numérique des trapèzes : (a-b)/n * (f(x_0)+2f(x_1)+2f(x_2)+...+f(x_n))
     * où a et b sont les bornes de l'intégrales et n est le nombre de points de discrétisations. 
     * Ici, (a+b)/n = pas de discrétisation selon x (=l) ou y (=h).
 * Arguments
 * =========
 * cx (input) - indices du contour selon x
 * cy (input) - indices du contour selon y
 * grid (input) - paramètres de la grille
 * physical (input) - paramètres physiques
 * output (input) - données de sortie
 */
{
    printf("\n------ compute_charge.c -> compute_flux() ------\n");

    int i, j;
    int nx = grid->nx;
    int ny = grid->ny;
    double h = physical->Le / (grid->my - 1);
    double l = physical->le / (grid->mx - 1);
    double eps_0 = 8.854e-12, Q;
    output->flux = 0;

    if (cx[0] < 0 || cx[1] > nx || cy[0] < 0 || cy[1] > ny)
    {
        printf("Erreur: les indices du contour dépassent le domaine\n");
        return 1;
    }

    printf("- Contour : cx = [%d, %d], cy = [%d, %d]\n", cx[0], cx[1], cy[0], cy[1]);
    for (i=cx[0]; i<=cx[1]; i++)
    {
        /* On parcourt les indices du contour selon x et on ajoute la composante y des vecteur au flux */
        output->flux -= fabs(output->EField[i][cy[0]][1] + output->EField[i+1][cy[0]][1])*h/2; /* On soustrait les vecteurs du champ électrique pour le sud du contour puisque le produit scalaire du vecteur normal et du champ est négatif (directions opposées) */
        output->flux -= fabs(output->EField[i][cy[1]][1] + output->EField[i][cy[1]][1])*h/2; /* On additionne les vecteurs du champ électrique pour le nord du contour puisque le produit scalaire du vecteur normal et du champ est négatif (directions opposées) */
    }

    for (j=cy[0]; j<=cy[1]; j+=7)
    {
        /* On parcourt les indices du contour selon y et on ajoute la composante x des vecteur au flux */
        output->flux -= fabs(output->EField[cx[0]][j][0] + output->EField[cx[0]][j+1][0])*l/2; /* Même raison que pour la boucle for précédente */
        output->flux -= fabs(output->EField[cx[1]][j][0] + output->EField[cx[1]][j+1][0])*l/2; /* Même raison que pour la boucle for précédente */
    }

    printf("- Flux = %f\n", output->flux);
    Q = output->flux * eps_0;
    printf("- Q = %.15e [Coulomb/m]\n", Q);
    printf("- C = %f [µF/m]\n", (fabs(Q)/(physical->ue - physical->ui))*1e6);
    return 0;
}


