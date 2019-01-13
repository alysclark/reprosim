
#include "exports.h"
#include "utils.h"

#include <string.h>

void export_1d_elem_field_c(int *ne_field, const char *EXELEMFILE, int *EXELEMFILE_LEN,
                            const char *group_name, int *group_name_len, const char *field_name, int *field_name_len );
void export_1d_elem_geometry_c(const char *EXELEMFILE, int *EXELEMFILE_LEN, const char *name, int *name_len);
void export_node_field_c(int *nj_field, const char *EXNODEFIELD, int *EXNODEFIELD_LEN,
                         const char *name, int *name_len, const char *field_name, int *field_name_len);
void export_terminal_perfusion_c(const char *EXNODEFILE, int *EXNODEFILE_LEN, const char *name, int *name_len);
void export_node_geometry_c(const char *EXNODEFILE, int *EXNODEFILE_LEN, const char *name, int *name_len);
void export_cell_location_c(const char *filename,int *FILENAME_LEN, int *cell_population);
void export_cell_exnode_c(const char *filename,int *FILENAME_LEN, int *cell_population);

void export_1d_elem_field(int ne_field, const char *EXELEMFILE, const char *group_name, const char *field_name )
{
  int filename_len = strlen(EXELEMFILE);
  int group_name_len = strlen(group_name);
  int field_name_len = strlen(field_name);

  export_1d_elem_field_c(&ne_field, EXELEMFILE, &filename_len, group_name, &group_name_len, field_name, &field_name_len);
}

void export_1d_elem_geometry(const char *EXELEMFILE, const char *name)
{
  int filename_len = strlen(EXELEMFILE);
  int name_len = strlen(name);

  export_1d_elem_geometry_c(EXELEMFILE, &filename_len, name, &name_len);
}

void export_node_field(int nj_field, const char *EXNODEFIELD, const char *name, const char *field_name)
{
  int filename_len = strlen(EXNODEFIELD);
  int name_len = strlen(name);
  int field_name_len = strlen(field_name);

  export_node_field_c(&nj_field, EXNODEFIELD, &filename_len, name, &name_len, field_name, &field_name_len);
}

void export_terminal_perfusion(const char *EXNODEFILE, const char *name)
{
  int filename_len = strlen(EXNODEFILE);
  int name_len = strlen(name);

  export_terminal_perfusion_c(EXNODEFILE, &filename_len, name, &name_len);
}

void export_node_geometry(const char *EXNODEFILE, const char *name)
{
  int filename_len = strlen(EXNODEFILE);
  int name_len = strlen(name);

  export_node_geometry_c(EXNODEFILE, &filename_len, name, &name_len);
}

void export_cell_location(const char *filename, int cell_population)
{
  int filename_len = strlen(filename);

  export_cell_location_c(filename, &filename_len, &cell_population);
}


void export_cell_exnode(const char *filename, int cell_population)
{
  int filename_len = strlen(filename);

  export_cell_exnode_c(filename, &filename_len, &cell_population);
}
