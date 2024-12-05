# PROJET

Pour exécuter le code avec les paramètres par défaut, il suffit d'exécuter les commandes ci-dessus : 
```bash
$ make
$ ./main
```
Le graphique généré est enregistré dans le fichier `potential_and_field.png`.

## `main.c`
Pour changer les paramètres de l'exécution, il suffit de modifier le fichier `main.c` et de recompiler le code avec la commande `make`.

### Expliquation des paramètres : 
- `int refinement` : le nombre de raffinements de la grille. Pour exécuter le code avec la grille de base, il faut intialiser la variable à `0`. Le nombre de raffinements maximal que j'ai réussi à calculer sur mon ordinateur est `8`.
- `int disk` : Permet de choisir si on veut une densité de charge nulle partout (sauf aux bords) en intialisant `disk` à `1` ou si on veut une densité de charge nulle partout sauf sur le disque mentionné dans les consignes en initialisant `disk` à 0.
- `double Q` : la charge du disque.
- `int nice_scale` : permet de choisir si on veut que l'échelle de couleur soit définie en fonction de la valeur maximale de la grille ou si les valeurs de celle-ci sont fixées entre les conditions aux bords et dans le trou. Il est utile d'initialiser`nice_scale` à `1` lorsque `disk` est intialisé à 1 puisque les valeurs êxtrèmes de la grilles ne seront plus les conditions aux bords et dans le trou. Dans les autres cas, il vaut mieux initialiser `nice_scale` à `0` afin d'économiser du temps de calcul.
- `int umfpack` : Lorsqu'elle est initialisée à `1` -> le code utilise le solveur `umfpack` tandis que quand elle est initialisé à `0`, le code utilise le solveur `agmg`.
- `int cx` et `int cy` : les indices de début et de fin du contour qui permet de calculer le flux du champ électrique. Ils sont définis plus tard dans le fichier `main.c`.

## `compute_density.c`
Ce fichier contient les fonctions qui permettent d'intialiser le tableau de densité (`initialize_density()`) (qui sert ensuite à calculer la matrice b) et de le remplir suivant certaines conditions (`compute_density_disk()`).

## `prob.c`
Ce fichier contient la fonction `prob()` qui remplit les matrices A et B en fonction de la densité de charge et des conditions aux bords.

## `compute_residue.c`
Ce fichier contient la fonction `compute_residue()` qui permet de calculer la norme du résidu de la solution.

## `compute_charge.c`
Ce fichier contient les fonctions : 
- `compute_EField_x` et `compute_EField_y` qui permettent de calculer les composantes du champ électrique pour un point de la grille.
- `compute_EField` qui utilise les 2 fonctions précédentes pour calculer le champ électrique en chaque point de la grille. 
- `compute_flux` qui calcule via la méthode des trapèzes le flux du champ électrique à travers un contour en Vm.

## `visualize.c`
Ce fichier contient les fonctions qui permettent de visualiser les résultats :
- `write_data()` écrit le potentiel dans un fichier.
- `write_guplot_script()` écrit le script gnuplot qui permet de générer le graphique.

## Solveur
- `umfpack` : https://users.encs.concordia.ca/~krzyzak/R%20Code-Communications%20in%20Statistics%20and%20Simulation%202014/Zubeh%F6r/SuiteSparse/UMFPACK/Doc/QuickStart.pdf
