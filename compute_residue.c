#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "utils.h"

int compute_residue(OutputData *output)
/* But 
* ===
    * Calculer le résidu de la solution x
    * Résidu de la forme : || Au - b ||_2 / || b ||_2
* Arguments
* =========
* output (input) - données de sortie
*/
{
    printf("\n------ compute_residue.c -> compute_residue() ------\n");
    /*FILE *f = fopen("residue.txt", "w");*/ // -> test file to see the values of the residue
    int n = output->n;

    int row_start, row_end, i, j, col;
    long double sum, residue, numerator_norm, denominator_norm;
    long double *y = malloc(n * sizeof(long double));
    if (y == NULL)
    {
        printf("Erreur lors de l'allocation de y\n");
        return 1;
    }

    residue = 0;
    numerator_norm = 0;
    denominator_norm = 0;
    for (i=0; i < n; i++)
    {
        sum = 0.0;
        row_start = (output->ia)[i]; // Le début de la ligne i est à l'indice ia[i]
        row_end = (output->ia)[i+1]; // La fin de la ligne i est à l'indice ia[i+1]

        for (j=row_start; j < row_end; j++)
        {
            col = (output->ja)[j]; // La colonne j de la ligne i est à l'indice ja[j]
            sum += (output->a)[j] * (output->x)[col]; // On ajoute le produit de la valeur de la matrice A et de la valeur de x (à la ligne i et a la colonne col) à la somme
        }
        y[i] = sum - output->b[i]; // Au - b 
        numerator_norm += y[i] * y[i]; // || Au - b ||_2
        denominator_norm += output->b[i] * output->b[i]; // || b ||_2
        /*fprintf(f, "sum = %Lf, b[%d] = %f, y[%d] = %Lf\n", sum, i, (output->b[i]), i,  y[i]);*/ // -> print the values of sum, b[i] and y[i] in the residue.txt file
    }

    if (denominator_norm == 0)
    {
        printf("Erreur : le vecteur b est nul\n");
        return 1;
    }
    residue = sqrt(numerator_norm) / sqrt(denominator_norm); // || Au - b ||_2 / || b ||_2
    printf("- residue = %.15Le\n", residue);

    free(y);
    return 0;
}
