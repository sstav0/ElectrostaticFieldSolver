#include "utils.h"

int compute_EField(GridParameters *grid, PhysicalParameters *physical, OutputData *output);
double compute_EField_x(double u1, double u2, double h);
double compute_EField_y(double u1, double u2, double l);
int compute_flux(int cx[2], int cy[2], GridParameters *grid, PhysicalParameters *physical, OutputData *output);
