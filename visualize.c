#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "visualize.h"
#include "prob.h"
#include "utils.h"


int write_data(GridParameters *grid, OutputData *output, PhysicalParameters *physical, char *filename)
/* But
* ===
    * Ecrire le potentiel u dans un fichier
* Arguments
* =========
* grid (input) - paramètres de la grille
* output (input) - données de sortie
* physical (input) - paramètres physiques
* filename (input) - nom du fichier de sortie
*/
{
    printf("\n------ visualize.c -> write_potential_to_file() ------\n");
    double Le = physical->Le;
    double le = physical->le;

    double h = Le/(grid->my-1); /* calcul du pas de discrétisation selon y*/
    double l = le/(grid->mx-1); /* calcul du pas de discrétisation selon x*/

    if (strcmp(filename + strlen(filename) - 4, ".dat") != 0)
    {
        fprintf(stderr, "Error: Output file must be a .dat file.\n");
        return 1;
    }

    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
        printf("Erreur lors de l'ouverture du fichier %s\n", filename);
        return 1;
    }

    /* Écriture des données dans le fichier */
    fprintf(fp, "# x y u\n");
    for (int i = 0; i < grid->mx; i++)
    {
        for (int j = 0; j < grid->my; j++)
        {
            if (output->grid_values[i] == NULL)
            {
                printf("Error: grid_values is NULL for j = %d, i = %d\n", j, i);
                return 1;
            }
            if (i > 0 && j > 0 && i <= grid->nx && j <= grid->ny) /* On n'écrit pas le potentiel pour les bords (i=0 et j=0) puisqu'il est nul */
                                                                /* [nx] et [ny] sont les derniers indices du tableau du champ */
            {
                if (output->EField[i][j] == NULL)
                {
                    printf("Error: EField is NULL for j = %d, i = %d\n", j, i);
                    return 1;
                }
                fprintf(fp, "%f %f %f %f %f\n", i*l, j*h, (output->grid_values)[i][j], output->EField[i][j][0], output->EField[i][j][1]);
            }
            else
                fprintf(fp, "%f %f %f\n", i * l, j * h, (output->grid_values)[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    return 0;
}


int write_gnuplot_script(char *filename, char *outputfile, PhysicalParameters *physical, OutputData *output)
/* But
* ===
    * Ecrire un script gnuplot pour visualiser le potentiel
* Arguments
* =========
* filename (input) - nom du fichier qui contient les données du potentiel²
*/
{
    printf("\n------ visualize.c -> write_gnuplot_script() ------\n");
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe == NULL)
    {
        fprintf(stderr, "Error opening pipe to Gnuplot.\n");
        return 1;
    }
    if (strcmp(outputfile + strlen(outputfile) - 4, ".png") != 0)
    {
        fprintf(stderr, "Error: Output file must be a .png file.\n");
        return 1;
    }

    fprintf(gnuplotPipe, "set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n"); // On 
    fprintf(gnuplotPipe, "set output '%s'\n", outputfile);
    fprintf(gnuplotPipe, "set title 'Electrostatic Potential [V] and Electrical Field [Vm]'\n");
    fprintf(gnuplotPipe, "set xlabel 'X-axis [mm]'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y-axis [mm]'\n");
    fprintf(gnuplotPipe, "unset key\n"); // Remove legend for cleaner visualization
    fprintf(gnuplotPipe, "set cblabel 'Potential [V]'\n");
    fprintf(gnuplotPipe, "set grid\n");
    fprintf(gnuplotPipe, "set tics nomirror\n");
    fprintf(gnuplotPipe, "set cbrange [%f:%f]\n", output->min, output->max);
    fprintf(gnuplotPipe, "set palette rgb 33,13,10\n");
    fprintf(gnuplotPipe, "set xrange [0:%f]\n", physical->le);
    fprintf(gnuplotPipe, "set yrange [0:%f]\n", physical->Le);

    fprintf(gnuplotPipe, "set bmargin at screen 0.15\n"); /* Décalage du graphique vers le haut */
    fprintf(gnuplotPipe, "set label 'Rq.: Les vecteurs du champ électrique sont normalisés et leur nouvelle norme est ensuite multipliée par 1/10' at screen 0.5, 0.04 center\n"); /* Ajout d'une note pour expliquer la normalisation des vecteurs */

    fprintf(gnuplotPipe, "plot '%s' with image, '%s' using 1:2:(($4/sqrt(($4)**2+($5)**2))*1e-1):($5/sqrt(($4)**2+($5)**2)*1e-1) every (2**(%d/1.5)):(2**(%d/1.5)) with vectors lc -1 filled \n", filename, filename, physical->refinement+1, physical->refinement+1);

    fflush(gnuplotPipe); // Ensure all commands are sent
    pclose(gnuplotPipe);
    return 0;
}
