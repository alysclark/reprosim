#include "abm_models.h"
#include "string.h"

void initialise_abm_c(const char *model_type,  int *model_type_len, int *total_cells, int *num_forces, double *time_step, double *min_time_step);
void move_cells_force_c(int *cell_population, double *kdrag, double *input_dt);
void check_cell_tube_c(int *cell_population);
double get_current_t_c();
double get_current_dt_c();

void initialise_abm(const char *model_type, int total_cells, int num_forces, double time_step, double min_time_step)
{
    int model_type_len = strlen(model_type);
    initialise_abm_c(model_type, &model_type_len, &total_cells, &num_forces, &time_step, &min_time_step);
}

void move_cells_force(int cell_population, double kdrag, double input_dt)
{
	move_cells_force_c(&cell_population, &kdrag, &input_dt);
}

void check_cell_tube(int cell_population)
{
	check_cell_tube_c(&cell_population);
}

double get_current_t()
{
  return get_current_t_c();
}

double get_current_dt()
{
  return get_current_dt_c();
}
