#include <stdio.h>
#include "print_matrix.h"
#include "utils.h"

void print_matrix(OutputData *output)
{
/*
 * But
   ===
   Afficher une matrice stockée dans le format CRS défini par les trois vecteurs ia, ja et a (où a est le vecteur qui contient les données)

   Arguments
   =========
   n (input) - nombre de lignes de la matrice
   nnz (input) - nombre de colonnes de la matrice
   ia    (input) - vecteur des débuts de lignes
   ja    (input) - vecteur des indices de colonnes
   a     (input) - vecteur des valeurs non nulles
*/
    int n = output->n;
    int nnz = output->nnz;

    int i, j, nz_index;
    FILE *f = fopen("matrix.txt", "w");
    FILE *fia = fopen("ia.txt", "w");
    FILE *fja = fopen("ja.txt", "w");
    FILE *fa = fopen("a.txt", "w");
    FILE *fb = fopen("b.txt", "w");
    FILE *fx = fopen("x.txt", "w");

    // Loop over each row
    for (i = 0; i < n; i++) {
        nz_index = (output->ia)[i];  // Start index of non-zero elements in this row

        // Loop over each column
        for (j = 0; j < nnz; j++) {
            // Check if current column has (output->a) non-zero element
            if (nz_index < (output->ia)[i+1] && (output->ja)[nz_index] == j) {
                fprintf(f, "%6.2f ", (output->a)[nz_index]);  // Print non-zero element
                /*printf("%6.2f ", (output->a)[nz_index]);  // Print non-zero element*/
                nz_index++;  // Move to the next non-zero element
            } else {
                fprintf(f, "%6.2f ", 0.0);  // Print zero for empty positions
                /*printf("%6.2f ", 0.0);  // Print zero for empty positions*/
            }
        }
        fprintf(f, "\n");  // Newline after each row
        /*printf("\n");  // Newline after each row*/
    }

    // Print the (output->ia) vector
    for (i = 0; i < n+1; i++) {
        fprintf(fia, "%d\n", (output->ia)[i]);
        /*printf("%d ", (output->ia)[i]);*/
    }

    // Print the (output->ja) vector
    for (i = 0; i < nnz; i++) {
        fprintf(fja, "%d\n", (output->ja)[i]);
        /*printf("%d ", (output->ja)[i]);*/
    }

    // Print the (output->a) vector
    for (i = 0; i < nnz; i++) {
        fprintf(fa, "%f\n", (output->a)[i]);
        /*printf("%6.2f ", (output->a)[i]);*/
    }

    // Print the (output->b) vector
    for (i = 0; i < n; i++) {
        fprintf(fb, "%f\n", (output->b)[i]);
        /*printf("%6.2f ", (output->b)[i]);*/
    }

    for (i = 0; i < n; i++) {
        fprintf(fx, "%f\n", (output->x)[i]);
        /*printf("%6.2f ", (output->x)[i]);*/
    }
    fclose(f); fclose(fia); fclose(fja); fclose(fa); fclose(fb); fclose(fx);
}
