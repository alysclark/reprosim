module mix_fem_c
implicit none
private

contains

!!!###################################################################################

subroutine read_b_matrix_c(filename,filename_len) bind(C, name="read_b_matrix_c")

use iso_c_binding, only: c_ptr
use utils_c, only: strncpy
use other_consts, only: MAX_FILENAME_LEN
use mix_fem, only: read_b_matrix
implicit none

type(c_ptr), value, intent(in) :: filename
integer,intent(in) :: filename_len
character(len=MAX_FILENAME_LEN) :: filename_f

call strncpy(filename_f, filename, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
call so_read_b_matrix(filename_f,filename_len)
#else
call read_b_matrix(filename_f,filename_len)
#endif

end subroutine read_b_matrix_c

!!!###################################################################################

subroutine read_e2face_c(filename,filename_len) bind(C, name="read_e2face_c")

use iso_c_binding, only: c_ptr
use utils_c, only: strncpy
use other_consts, only: MAX_FILENAME_LEN
use mix_fem, only: read_e2face
implicit none

type(c_ptr), value, intent(in) :: filename
integer,intent(in) :: filename_len
character(len=MAX_FILENAME_LEN) :: filename_f

call strncpy(filename_f, filename, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
call so_read_e2face(filename_f,filename_len)
#else
call read_e2face(filename_f,filename_len)
#endif

end subroutine read_e2face_c


!!!###################################################################################

subroutine read_face2e_c(filename,filename_len) bind(C, name="read_face2e_c")

use iso_c_binding, only: c_ptr
use utils_c, only: strncpy
use other_consts, only: MAX_FILENAME_LEN
use mix_fem, only: read_face2e
implicit none

type(c_ptr), value, intent(in) :: filename
integer,intent(in) :: filename_len
character(len=MAX_FILENAME_LEN) :: filename_f

call strncpy(filename_f, filename, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
call so_read_face2e(filename_f,filename_len)
#else
call read_face2e(filename_f,filename_len)
#endif

end subroutine read_face2e_c

!!!###################################################################################

subroutine assemble_sparse_matrices_c() bind(C, name="assemble_sparse_matrices_c")

use mix_fem, only: assemble_sparse_matrices
implicit none


#if defined _WIN32 && defined __INTEL_COMPILER
call so_assemble_sparse_matrices()
#else
call assemble_sparse_matrices()
#endif

end subroutine assemble_sparse_matrices_c

!!!###################################################################################

subroutine create_sampling_grid_c() bind(C, name="create_sampling_grid_c")

use mix_fem, only: create_sampling_grid
implicit none


#if defined _WIN32 && defined __INTEL_COMPILER
call so_create_sampling_grid()
#else
call create_sampling_grid()
#endif

end subroutine create_sampling_grid_c


end module mix_fem_c
