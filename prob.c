#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "print_matrix.h"
#include "compute_density.h"
#include "utils.h"
#include "visualize.h"
#include "prob.h"

int prob(double **rho, GridParameters *grid, PhysicalParameters *physical, OutputData *output)
/*
   But
   ===
   Générer le système linéaire n x n 
                          Au = b
   qui correspond à la discrétisation sur une grille cartésienne 
   régulière m x m de l'équation de Poisson à deux dimensions

            d    d        d    d
         - == ( == u ) - == ( == u )  = 0     sur [0,1] x [0,1]
           dx   dx       dy   dy

  avec les conditions aux limites de Dirichlet
         u = 1  sur (x,0), (x,1), (0,y) et (1,y), avec 0 <= x,y <= 1 . 

  La numérotation des inconnues est lexicographique, la direction x étant 
  parcourue avant celle de y. La matrice A est retournée dans le format 
  CRS qui est défini via trois tableaux : 'ia', 'ja' et 'a'.

  N.B. Ce problème possède une solution analytique assez simple.

  Arguments
  =========
  rho (input) - tableau de taille n contenant la densité de charge pour chaque variable du domaine

  grid(input) (voir fichier "utils.h") : 
    - mx (input)  - nombre de points dans la direction x de la grille (mw <= 2 n'est pas valide)
    - my (input)  - nombre de points dans la direction y de la grille (my <= 2 n'est pas valide)
    - nx (output) - pointeur vers le nombre d'inconnues dans la direction x
    - ny (output) - pointeur vers le nombre d'inconnues dans la direction y

  physical(input) (voir fichier "utils.h") :
    - fac_eps (i-nput) - facteur de epsilon_0
    - ue (input)  - valeur de la solution analytique aux bords du domaine
    - ui (input)  - valeur de la solution aux bords du trou
    - Le (input)  - Largeur du domaine (Le = 0 n'est pas valide) 
    - le (input)  - intervalle en x (= longueur) du trou (double le[2] = {0,0} correspond à un trou inexistant si double le[2] = {0,0})
    - li (input)  - intervalle en x (= longueur) du trou (double li[2] = {0,0} correspond à un trou inexistant si double Li[2] = {0,0})
    - Li (input)  - intervalle en y (= largeur) du trou (double Li[2] = {0,0} correspond à un trou inexistant si double Li[2] = {0,0})

  output(output) (voir fichier "utils.h") :
    - n  (output) - Nombre de points de la grille
    - i_start_x (output) - Indice du point qui se situe en bas à gauche du trou
    - i_end_x (output) - Indice du point qui se situe en bas à droite du trou
    - i_start_y (output) - Indice du point qui se situe en bas à gauche du trou
    - i_end_y (output) - Indice du point qui se situe en haut à droite du trou
    - ia (output) - Tableau des indices de début de ligne
    - ja (output) - Tableau des indices de colonne
    - a (output) - Tableau des valeurs non nulles de la matrice A
    - b (output) - Tableau des valeurs de la matrice B
    - x (output) - Tableau des valeurs de la matrice de solutions
    - grid_values (output) - Tableau des valeurs pour chaque point de la grille
    - EField (output) - Tableau des valeurs du champ électrique pour chaque point de la grille
    - flux (output) - Flux du champ électrique à travers un contour englobant le trou
    - min (output) - Valeur minimale du potentiel électrique de la grille
    - max (output) - Valeur maximale du potentiel électrique de la grille
*/
{
    printf("\n------ prob.c -> prob() ------\n");
    int ix, iy, ind, n_hole, i_start_x, i_end_x, i_start_y, i_end_y, nnz;
    int nx_hole = 0;
    int ny_hole = 0;
    double invh2, h, l, eps;

    int mx = grid->mx;
    int my = grid->my;
    int nx = grid->nx;
    int ny = grid->ny;

    int ue = physical->ue;
    int ui = physical->ui;
    double le = physical->le;
    double Le = physical->Le;
    double li[2] = {physical->li[0], physical->li[1]};
    double Li[2] = {physical->Li[0], physical->Li[1]};
    double fac_eps = physical->fac_eps;

    h = Le/(my-1); /* calcul du pas de discrétisation selon y*/
    l = le/(mx-1); /* calcul du pas de discrétisation selon x*/

    printf("- Pas de discrétisations : h = %f [mm], l = %f [mm]\n", h, l);
    eps = fac_eps * 8.854e-12; /* epsilon_0 */

    if(nx < 0) {
        printf("\n ERREUR : mx = %d n'est pas une valeur valide\n\n",mx);
        return 1;
    }
    if (ny < 0) {
        printf("\n ERREUR : mx = %d n'est pas une valeur valide\n\n",my);
        return 1;
    }

    invh2 = 1/(h*h) * 1e6; /* h^-2 pour L != 1 -> pourrait être simplifié */

    /* Calcul du nombre d'inconnues en fonction des dimensions du trou */

    output->n  = nx * ny; /* nombre d'inconnues sans compter le trou */
    printf("- Total de points n = %d\n", output->n);

    if ((li[0] == 0 && li[1] == 0) && (Li[0] == 0 && Li[1] == 0)) {
        /* No hole is defined, so output->n remains unchanged */
        n_hole = 0;
        i_start_x = 0;
        i_start_y = 0;
        i_end_x = 0;
        i_end_y = 0; 

    } else if ((li[0] == 0 && li[1] == 0) != (Li[0] == 0 && Li[1] == 0)) {
        printf("\n ERREUR : le trou est mal défini\n\n");
        return 1;
    }else /* il y a un trou et il est bien défini */
    {
        /* Calcul du nbre d'inconnues en fonctions des dimensions du trou */
        i_start_x = ceil(li[0] / h) - 1; /* ceil pour arrondir à l'entier supérieur et - 1 pour que les indices aient 0,0 comme repère*/
        i_end_x = floor(li[1] / h) - 1; /* floor pour arrondir à l'entier inférieur */
        i_start_y = ceil(Li[0] / l) - 1;
        i_end_y = floor(Li[1] / l) - 1;

        printf("- Indices du trou : i_start_x = %d, i_end_x = %d, i_start_y = %d, i_end_y = %d\n", i_start_x, i_end_x, i_start_y, i_end_y);

        nx_hole = i_end_x - i_start_x + 1;
        ny_hole = i_end_y - i_start_y + 1;

        n_hole = nx_hole * ny_hole;

        output->n = output->n - n_hole; /* nombre d'inconnues en enlevant celles dans le trou */
        printf("- Nombre de points dans le trou : %d\n", n_hole);
    }

    /* Remplissage du tableau vertices_hole pour les autres fonctions du fichier main.c*/
    output->i_start_x = i_start_x;
    output->i_end_x = i_end_x;
    output->i_start_y = i_start_y;
    output->i_end_y = i_end_y;

    /* Je passe un pointeur vers nnz à compute_nnz pour modifier la variable nnz directement dans la fonction */
    if (compute_nnz(n_hole, h, grid, output, physical))
    {
        printf("\n ERREUR : problème dans le calcul du nombre d'éléments non nuls de la matrice A\n\n");
        return 1;
    }
    printf("- Nombre d'éléments non nuls de la matrice A (fonction compute_nnz) : %d\n", output->nnz);

    /* allocation des tableaux */

    output->ia  = malloc((output->n + 1) * sizeof(int));
    output->ja  = malloc(((output->nnz)) * sizeof(int));
    output->a   = malloc(((output->nnz)) * sizeof(double));
    output->b   = malloc(output->n * sizeof(double));

    /* allocation réussite? */

    if (output->ia == NULL || output->ja == NULL || output->a == NULL || output->b == NULL ) {
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n\n");
        return 1;
    }

    /* partie principale : remplissage de la matrice */
    ind = 0;
    nnz = 0;
    n_hole = 0;
    for (iy = 0; iy < ny; iy++) {
        for (ix = 0; ix < nx; ix++) {

            if ((ix >= i_start_x && ix <= i_end_x) && (iy >= i_start_y && iy <= i_end_y))
            {
                n_hole += 1; /* on compte le nombre de points dans le trou */
                continue;
            }

            /* numéro de l'équation (l'indice 0 représente la 1e inconnue en bas à gauche) */
            ind = ix + nx * iy - n_hole;

            /* marquer le début de la ligne suivante dans le tableau 'ia' */
            (output->ia)[ind] = nnz;

            /* calculer le membre de droite */
            (output->b)[ind] = rho[ix][iy]/eps; /* rho contient la densité de chaque point de la grille */

            /* remplissage de la ligne : voisin sud */
            if (iy > 0)
            {
                /* Si on se trouve dessu du trou alors le voisin sud est une condition au bord du trou */
                if (iy == i_end_y + 1 && ix >= i_start_x && ix <= i_end_x) 
                    (output->b)[ind] += ui * invh2;
                else
                {
                    (output->a)[nnz] = -invh2;
                    /* Calcul de l'indice du voisin sud */
                    if (iy >= i_start_y + 1 && iy <= i_end_y) /* -> Si on se trouve à la même hauteur (y) que le trou */
                        (output->ja)[nnz] = ind - nx + nx_hole;
                    else if (iy == i_start_y && ix >= i_end_x) /* -> Si on se trouve à droite du trou et au niveau de sa 1e ligne */
                        (output->ja)[nnz] = ind - nx + nx_hole;
                    else if (iy == i_end_y + 1 && ix<= i_end_x) /* -> Si on se trouve partout sauf à droite du trou et au dessus de sa dernière ligne */
                        (output->ja)[nnz] = ind - nx + nx_hole;
                    else
                        (output->ja)[nnz] = ind - nx;
                    nnz++;
                }
            }
            else
                (output->b)[ind] += ue * invh2; /* condition de Dirichlet à la frontière sud */

            /* remplissage de la ligne : voisin ouest */
            if (ix > 0)
            {
                /* Si on se trouve à droite du trou alors le voisin ouest est une condition au bord du trou */
                if (ix == i_end_x + 1 && iy >= i_start_y && iy <= i_end_y) // 
                    (output->b)[ind] += ui * invh2;
                else
                {
                    (output->a)[nnz] = -invh2;
                    (output->ja)[nnz] = ind - 1;
                    nnz++;
                }
            }
            else
                (output->b)[ind] += ue * invh2; /* condition de Dirichlet à la frontière ouest */

            /* remplissage de la ligne : élément diagonal */
            (output->a)[nnz] = 4.0*invh2; 
            (output->ja)[nnz] = ind;
            nnz++;

            /* remplissage de la ligne : voisin est */
            if (ix < nx - 1)
            {
                if (ix == i_start_x - 1 && iy >= i_start_y && iy <= i_end_y)
                    (output->b)[ind] += ui * invh2;
                else
                {
                    (output->a)[nnz] = -invh2; 
                    (output->ja)[nnz] = ind + 1;
                    nnz++;
                }
            }
            else
                (output->b)[ind] += ue * invh2; /* condition de Dirichlet à la frontière est */

            /* remplissage de la ligne : voisin nord */
            if (iy < ny - 1)
            {
                /* Si on se trouve en dessous du trou alors le voisin nord est une condition au bord du trou */
                if (iy == i_start_y - 1 && ix >= i_start_x && ix <= i_end_x)
                    (output->b)[ind] += ui * invh2;
                else
                {
                    (output->a)[nnz] = -invh2;
                    if (iy >= i_start_y && iy <= i_end_y -1)
                        (output->ja)[nnz] = ind + nx - nx_hole;
                    else if (iy == i_start_y - 1 && ix > i_end_x)
                        (output->ja)[nnz] = ind + nx - nx_hole;
                    else if (iy == i_end_y && ix < i_start_x)
                        (output->ja)[nnz] = ind + nx - nx_hole;
                    else
                        (output->ja)[nnz] = ind + nx;
                    nnz++;
                }
            }
            else
                (output->b)[ind] += ue * invh2; // condition de Dirichlet à la frontière nord
        }
    }

    /* dernier élément du tableau 'ia' */
    (output->ia)[ind + 1] = nnz;

    if (nnz != output->nnz)
    {
        printf("\n ERREUR : nnz et output->nnz ne correspondent pas. Le nombre d'éléments non nuls de la matrice A est incorrect\n\n");
        return 1;
    }

    return 0;
}
