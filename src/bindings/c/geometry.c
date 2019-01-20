
#include "geometry.h"

#include "string.h"

void add_matching_mesh_c();
void append_units_c();
void calc_capillary_unit_length_c(int *num_convolutes, int *num_generations);
void define_1d_elements_c(const char *ELEMFILE, int *filename_len);
void define_node_geometry_c(const char *NODEFILE, int *filename_len);
void define_rad_from_geom_c(const char *order_system, int *order_system_len, double *control_param,
                            const char *start_from, int *start_from_len, double *start_rad,
                            const char *group_type, int *group_type_len, const char *group_options, int *group_options_len);
void element_connectivity_1d_c();
void evaluate_ordering_c();
void read_icem_msh_c(const char *filename, int *filename_len);
void read_k_file_c(const char *filename, int *filename_len);

void add_matching_mesh()
{
	add_matching_mesh_c();
}

void append_units()
{
  append_units_c();
}

void calc_capillary_unit_length(int num_convolutes, int num_generations)
{
  calc_capillary_unit_length_c(&num_convolutes,&num_generations);
}

void define_1d_elements(const char *ELEMFILE)
{
  int filename_len = strlen(ELEMFILE);
  define_1d_elements_c(ELEMFILE, &filename_len);
}

void define_node_geometry(const char *NODEFILE)
{
  int filename_len = strlen(NODEFILE);
  define_node_geometry_c(NODEFILE, &filename_len);
}

void define_rad_from_geom(const char *order_system, double control_param, const char *start_from,
                          double start_rad, const char*group_type, const char *group_options)
{
  int order_system_len = strlen(order_system);
  int start_from_len = strlen(start_from);
  int group_type_len = strlen(group_type);
  int group_options_len = strlen(group_options);
  define_rad_from_geom_c(order_system, &order_system_len, &control_param, start_from, &start_from_len, &start_rad,
                         group_type, &group_type_len, group_options, &group_options_len);

}

void element_connectivity_1d()
{
  element_connectivity_1d_c();
}

void evaluate_ordering()
{
  evaluate_ordering_c();
}

void read_icem_msh(const char *filename)
{
  int filename_len = strlen(filename);
  read_icem_msh_c(filename, &filename_len);
}

void read_k_file(const char *filename)
{
  int filename_len = strlen(filename);
  read_k_file_c(filename, &filename_len);
}
