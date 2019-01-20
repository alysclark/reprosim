
#ifndef REPROSIM_GEOMETRY_H
#define REPROSIM_GEOMETRY_H

#include "symbol_export.h"

SHO_PUBLIC void add_matching_mesh();
SHO_PUBLIC void append_units();
SHO_PUBLIC void calc_capillary_unit_length(int num_convolutes, int num_generations);
SHO_PUBLIC void define_1d_elements(const char *ELEMFILE);
SHO_PUBLIC void define_node_geometry(const char *NODEFILE);
SHO_PUBLIC void define_rad_from_geom(const char *ORDER_SYSTEM, double CONTROL_PARAM, const char *START_FROM,
                                     double START_RAD, const char *GROUP_TYPE, const char *GROUP_OPTIONS);
SHO_PUBLIC void element_connectivity_1d();
SHO_PUBLIC void evaluate_ordering();
SHO_PUBLIC void read_icem_msh(const char *filename);
SHO_PUBLIC void read_k_file(const char *filename);

#endif /* REPROSIM_GEOMETRY_H */
