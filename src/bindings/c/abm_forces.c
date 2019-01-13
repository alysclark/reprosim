#include "abm_forces.h"
#include "string.h"

void calc_random_forces_c(int *cell_population,  int *force_field, double *force_magnitude);
void calc_saghian_chemo_forces_c(int *cell_population,  int *force_field, double *fradial, double *cradial, double *faxial, double *caxial);

void calc_random_forces(int cell_population,  int force_field, double force_magnitude)
{
    calc_random_forces_c(&cell_population, &force_field, &force_magnitude);
}

void calc_saghian_chemo_forces(int cell_population,  int force_field, double fradial, double cradial, double faxial, double caxial)
{
	calc_saghian_chemo_forces_c(&cell_population, &force_field, &fradial, &cradial, &faxial, &caxial);
}
