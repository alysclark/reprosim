#ifndef AETHER_MIX_FEM_H
#define AETHER_MIX_FEM_H

#include "symbol_export.h"

SHO_PUBLIC void read_b_matrix(const char *filename);
SHO_PUBLIC void read_e2face(const char *filename);
SHO_PUBLIC void read_face2e(const char *filename);
SHO_PUBLIC void assemble_sparse_matrices();
SHO_PUBLIC void create_sampling_grid();
SHO_PUBLIC void resample_grid();
SHO_PUBLIC void compute_body_forces(double inletPressure, double outletPressure);
SHO_PUBLIC void define_velocity_at_cell(int ccount, double velocity_at_cell, int ijk);

#endif /* AETHER_MIX_FEM_H */
