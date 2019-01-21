#include "mix_fem.h"
#include "string.h"

void read_b_matrix_c(const char *filename,  int *filename_len);
void read_e2face_c(const char *filename,  int *filename_len);
void assemble_sparse_matrices_c();
void create_sampling_grid_c();

void read_b_matrix(const char *filename)
{
    int filename_len = strlen(filename);
    read_b_matrix_c(filename, &filename_len);
}

void read_e2face(const char *filename)
{
    int filename_len = strlen(filename);
    read_e2face_c(filename, &filename_len);
}


void assemble_sparse_matrices()
{
    assemble_sparse_matrices_c();
}

void create_sampling_grid()
{
	create_sampling_grid_c();
}
