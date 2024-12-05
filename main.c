#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "prob.h"
#include "umfpk.h"
#include "time.h"
#include "print_matrix.h"
#include "prob.h"
#include "compute_density.h"
#include "compute_residue.h"
#include "visualize.h"
#include "utils.h"
#include "compute_charge.h"

/* Déclaration de la fonction dagmg_ */
void dagmg_(int*,double*,int*,int*,double*,double*,int*,int*,int*,int*,double*);

/* Fonction main */
int main(int argc, char *argv[])
{
    /* ------------ VARIABLES A MODIFIER ------------ */
    int refinement = 2; /* raffinement de la grille -> 0 correspond à la grille de base */
    int disk = 1; /* 1 pour un disque, 0 pour rho = 0 partout */
    double Q = 1e-2; /* Charge du disque */
    int nice_scale = 1; /* 1 pour une échelle de couleur en fonction de la valeur maximale de la grille (utile quand la densité de charge n'est pas nulle partout dans les conditions initiales) et 0 pour que les valeurs max et min de l'échelle soient les conditions initiales de potentiel */

    /* ------------ VARIABLES A NE PAS MODIFIER (sauf pour changer la physique du problème) ------------ */
    /* déclarer les variables */
    int i, j;
    /* relatives à la grille */
    int mx = 10, my=9;
    int nx, ny;
    double le = 4.5, Le =4;
    double li[2] = {1, 3.5}; /* coordonnées de la longueur selon x du trou */
    double Li[2] = {1, 2}; /* coordonnées de la longueur selon y du trou */

    /*relative à la physique du problème */
    int ue = 5, ui=-5;
    double fac_eps = 8.5;

    /* relatives à la solution */
    double **rho;
    double tc1, tc2, tw1, tw2;

    /* relatives au contour */
    int cx[2], cy[2];

    /*Variables utiles à la visualisation */
    char *filename_data = "potential.dat";
    char *outputfile = "potential_and_field.png";

    /* Application du rafinement */
    mx = mx * pow(2,refinement) - (pow(2, refinement)-1); /* raffinement de la grille */
    my = my * pow(2,refinement) - (pow(2, refinement)-1); /* raffinement de la grille */
    nx = mx - 2;
    ny = my - 2;

    /* Initialisation des structures */
    GridParameters grid = {mx, my, nx, ny};
    PhysicalParameters physical = {ue, ui, le, Le, {li[0], li[1]}, {Li[0], Li[1]},  fac_eps, refinement}; 
    OutputData output;

    /* Allouer de la mémoire pour le vecteur de densité */
    rho = (double **)malloc(nx * sizeof(double *));
    for (int i = 0; i < nx; i++)
        rho[i] = (double *)malloc(ny * sizeof(double));

    /* Calcul de la densité point par point */
    if (initialize_density(&rho, &grid))
        return 1;
    if (disk && compute_density_disk(&rho, Q, &grid, &physical)) /* Si disk est activé, on calcule la densité dans un disque */
        return 1;

    if (prob(rho, &grid, &physical, &output))
        return 1;
    printf("\n------ main.c -> main() ------");
    printf("\nPROBLEM:\n");
    printf("- Grid : mx = %d, my = %d\n- Physical : \n    ue (tension aux bord du domaine) = %d [V], ui (tension dans le trou) = %d [V]\n    le (longueur du domaine) = %f [mm], Le (largeur du domaine) = %f [mm]\n    fac_eps (facteur permittivité) = %f\n", grid.mx, grid.my, physical.ue, physical.ui, physical.le, physical.Le, physical.fac_eps);
    printf("- Output : \n    i_start_x = %d, i_end_x = %d, i_start_y = %d, i_end_y = %d\n    n (nbre d'inconnues) = %d, nx = %d, ny = %d\n", output.i_start_x, output.i_end_x, output.i_start_y, output.i_end_y, output.n, grid.nx, grid.ny);

    /* Libérer la mémoire inutile */
    for (i = 0; i < grid.nx; i++)
        free(rho[i]);
    free(rho);

    /* Allocation de mémoire pour le vecteur solution */
    output.x = malloc(output.n * sizeof(double));
    if (output.x == NULL ) {
        printf("\n ERREUR : pas de mémoire pour vecteur des solutions\n\n");
        return 1;
    }

    /* résoudre et mesurer le temps de solution */
    printf("\nSolving the system... -> solve_umfpack()");
    tc1 = mytimer_cpu(); tw1 = mytimer_wall();
    if(solve_umfpack(output.n, output.ia, output.ja, output.a, output.b, output.x))
        return 1;
    tc1 = mytimer_cpu(); tw1 = mytimer_wall();
    tc2 = mytimer_cpu(); tw2 = mytimer_wall();
    printf("\nTemps de solution (CPU): %5.1f sec",tc2-tc1);
    printf("\nTemps de solution (horloge): %5.1f sec \n",tw2-tw1);
    /*print_matrix(&output);*/

    /* Calcul du résidu */
    if (compute_residue(&output))
        return 1;

    /* Liberer la memoire */
    free(output.ia); 
    free(output.ja); 
    free(output.a); 
    free(output.b); 

    /* Calcul du champ électrique */
    if (compute_EField(&grid, &physical, &output))
        return 1;

    /* Choix du contour pour le calcul du flux */
    cx[0] = output.i_start_x - grid.nx/6;
    cx[1] = output.i_end_x + grid.nx/6;
    cy[0] = output.i_start_y - grid.nx/7;
    cy[1] = output.i_end_y + grid.nx/4;
    if (compute_flux(cx, cy, &grid, &physical, &output))
        return 1;

    /* Liberer la memoire inutile */
    free(output.x);

    /* Ecrire les données pour chaque point de la grille dans un fichier .dat */
    if (write_data(&grid, &output, &physical, filename_data))
        return 1;

    if (nice_scale)
        find_extremes_grid_values(&output, &grid);
    else
    {
        output.min = physical.ui;
        output.max = physical.ue;
    }

    /* Envoyer les commandes à gnuplot */
    if (write_gnuplot_script(filename_data, outputfile, &physical, &output))
        return 1;

    for (i = 0; i < mx; i++)
        free(output.grid_values[i]);
    free(output.grid_values);

    for (i=0; i <=nx; i++)
    {
        for (j=0; j <=ny; j++)
            free(output.EField[i][j]);
        free(output.EField[i]);
    }
    free(output.EField);
}
