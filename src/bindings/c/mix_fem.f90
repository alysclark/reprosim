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

subroutine resample_grid_c() bind(C, name="resample_grid_c")

use mix_fem, only: resample_grid
implicit none


#if defined _WIN32 && defined __INTEL_COMPILER
call so_resample_grid()
#else
call resample_grid()
#endif

end subroutine resample_grid_c

subroutine compute_body_forces_c(inletPressure,outletPressure) bind(C, name="compute_body_forces_c")
use arrays, only: dp
use mix_fem, only: compute_body_forces
implicit none
real(dp), intent(in):: inletPressure,outletPressure


#if defined _WIN32 && defined __INTEL_COMPILER
call so_compute_body_forces(inletPressure,outletPressure)
#else
call compute_body_forces(inletPressure,outletPressure)
#endif

end subroutine compute_body_forces_c



subroutine define_velocity_at_cell_c(ccount,velocity_at_cell,ijk) bind(C, name="define_velocity_at_cell_c")
use arrays, only: dp
use mix_fem, only: define_velocity_at_cell
implicit none
integer, intent(in) :: ccount,ijk
real(dp), intent(in):: velocity_at_cell


#if defined _WIN32 && defined __INTEL_COMPILER
call so_define_velocity_at_cell(ccount,velocity_at_cell,ijk)
#else
call define_velocity_at_cell(ccount,velocity_at_cell,ijk)
#endif

end subroutine define_velocity_at_cell_c


end module mix_fem_c
