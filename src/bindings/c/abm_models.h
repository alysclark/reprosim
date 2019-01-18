#ifndef AETHER_ABM_MODELS_H
#define AETHER_ABM_MODELS_H

#include "symbol_export.h"

SHO_PUBLIC void initialise_abm(const char *model_type, int total_cells, int num_forces, double time_step, double min_time_step);
SHO_PUBLIC void move_cells_force(int cell_population, double kdrag, double input_dt);
SHO_PUBLIC void check_cell_tube(int cell_population);

SHO_PUBLIC double get_current_t();
SHO_PUBLIC double get_current_dt();


#endif /* AETHER_ABM_MODELS_H */
