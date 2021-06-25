#ifndef REPROSIM_PRESSURE_RESISTANCE_FLOW_H
#define REPROSIM_PRESSURE_RESISTANCE_FLOW_H

#include "symbol_export.h"

SHO_PUBLIC void calculate_stats(const char *FLOW_GEN_FILE, double image_voxel_size);
SHO_PUBLIC void evaluate_prq(const char *mesh_type, const char *bc_type, const char *rheology_type, const char *vessel_type, double inlet_flow, double inlet_pressure, double outlet_pressure, double occlude_frac, int occlude_order, int occlude_side);


#endif /* REPROSIM_PRESSURE_RESISTANCE_FLOW_H */
