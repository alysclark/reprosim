#ifndef AETHER_ABM_FORCES_H
#define AETHER_ABM_FORCES_H

#include "symbol_export.h"

SHO_PUBLIC void calc_random_forces(int cell_population, int force_field, double force_magnitude);
SHO_PUBLIC void calc_saghian_chemo_forces(int cell_population, int force_field, double fradial, double cradial, double faxial, double caxial);
SHO_PUBLIC void calc_saghian_cell_cell(int cell_population, int force_field, double r0, double r1, double a, double b);

#endif /* AETHER_ABM_FORCES_H */
