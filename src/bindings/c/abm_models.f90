module abm_models_c
implicit none
private

contains

!!!###################################################################################

subroutine initialise_abm_c(model_type,model_type_len,total_cells,num_forces, &
    time_step,min_time_step) bind(C, name="initialise_abm_c")

use iso_c_binding, only: c_ptr
use utils_c, only: strncpy
use other_consts, only: MAX_STRING_LEN
use arrays, only: dp
use abm_models, only: initialise_abm
implicit none

type(c_ptr), value, intent(in) :: model_type
integer,intent(in) :: model_type_len, total_cells,num_forces
real(dp), intent(in) :: time_step, min_time_step
character(len=MAX_STRING_LEN) :: model_type_f

call strncpy(model_type_f, model_type, model_type_len)

#if defined _WIN32 && defined __INTEL_COMPILER
call so_initialise_abm(model_type_f,total_cells,num_forces,time_step,min_time_step)
#else
call initialise_abm(model_type_f,total_cells,num_forces,time_step,min_time_step)
#endif

end subroutine initialise_abm_c

!###################################################################################

subroutine move_cells_force_c(cell_population,kdrag,input_dt) bind(C, name = "move_cells_force_c")

use arrays, only:dp
use abm_models, only: move_cells_force

integer, intent(in) :: cell_population
real(dp), intent(in) :: kdrag
real(dp), intent(in) :: input_dt

#if defined _WIN32 && defined __INTEL_COMPILER
call so_move_cells_force_c(cell_population,kdrag,input_dt)
#else
call move_cells_force(cell_population,kdrag,input_dt)
#endif

end subroutine move_cells_force_c

!###################################################################################

subroutine check_cell_tube_c(cell_population) bind(C, name = "check_cell_tube_c")

use arrays, only:dp
use abm_models, only: check_cell_tube

integer, intent(in) :: cell_population

#if defined _WIN32 && defined __INTEL_COMPILER
call so_check_cell_tube_c(cell_population)
#else
call check_cell_tube(cell_population)
#endif

end subroutine check_cell_tube_c

!###################################################################################


function get_current_t_c() result(time_params) bind(C, name="get_current_t_c")
use arrays, only:dp
use abm_models, only: get_current_t
implicit none
real(dp) :: time_params
#if defined _WIN32 && defined __INTEL_COMPILER
time_params = so_get_current_t_c()
#else
time_params = get_current_t()
#endif
end function get_current_t_c

function get_current_dt_c() result(time_params) bind(C, name="get_current_dt_c")
use arrays, only:dp
use abm_models, only: get_current_dt
implicit none
real(dp) :: time_params
#if defined _WIN32 && defined __INTEL_COMPILER
time_params = so_get_current_dt_c()
#else
time_params = get_current_dt()
#endif
end function get_current_dt_c

end module abm_models_c
