module abm_forces_c
implicit none
private

contains

!!!###################################################################################

subroutine calc_random_forces_c(cell_population, force_field, force_magnitude) bind(C, name="calc_random_forces_c")
use arrays, only: dp
use abm_forces, only: calc_random_forces
implicit none

integer,intent(in) :: cell_population
integer, intent(in) :: force_field
real(dp), intent(in) :: force_magnitude

#if defined _WIN32 && defined __INTEL_COMPILER
call so_calc_random_forces(cell_population, force_field, force_magnitude)
#else
call calc_random_forces(cell_population, force_field, force_magnitude)
#endif

end subroutine calc_random_forces_c

!###################################################################################

  subroutine calc_saghian_chemo_forces_c(cell_population, force_field, fradial, cradial, faxial, caxial) &
       bind(C, name="calc_saghian_chemo_forces_c")
    use arrays, only: dp
    use abm_forces, only: calc_saghian_chemo_forces
    implicit none

    integer,intent(in) :: cell_population
    integer, intent(in) :: force_field
    real(dp), intent(in) :: fradial,cradial,faxial,caxial

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_calc_saghian_chemo_forces(cell_population, force_field,fradial,cradial,faxial,caxial)
#else
    call calc_saghian_chemo_forces(cell_population, force_field,fradial,cradial,faxial,caxial)
#endif

  end subroutine calc_saghian_chemo_forces_c

end module abm_forces_c
