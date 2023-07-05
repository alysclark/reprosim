
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

end module fetal_c