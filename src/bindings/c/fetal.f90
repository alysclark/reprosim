
module fetal_c
  implicit none
!Interfaces
private

contains
  
!
!> Perfusion fetal
  subroutine fetal_model_c() bind(C, name="fetal_model_c")

    use fetal, only: fetal_model
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_fetal_model()
#else
    call fetal_model()
#endif

  end subroutine fetal_model_c

  !
!> assign_fetal_arrays
  subroutine assign_fetal_arrays_c() bind(C, name="assign_fetal_arrays_c")

    use fetal, only: assign_fetal_arrays
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_assign_fetal_arrays()
#else
    call assign_fetal_arrays()
#endif

  end subroutine assign_fetal_arrays_c

end module fetal_c