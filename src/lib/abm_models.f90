module abm_models

  use solve, only: BICGSTAB_LinSolv,pmgmres_ilu_cr
  implicit none
  !Module parameters

  !Module types

  !Module depvar

  !Interfaces
  private
  public evaluate_abm
contains
!
!###################################################################################
!
subroutine evaluate_abm(model_type)
    use diagnostics, only: enter_exit,get_diagnostics_level
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_ABM" :: EVALUATE_ABM

    character(len=60), intent(in) :: model_type

    !local variables
    character(len=60) :: sub_name
    sub_name = 'evaluate_abm'
    call enter_exit(sub_name,1)
    write(*,*) model_type
    call allocate_abm_memory(model_type)
    call enter_exit(sub_name,2)
end subroutine evaluate_abm
!
!###################################################################################
!

subroutine allocate_abm_memory(model_type)
    use diagnostics, only: enter_exit,get_diagnostics_level
    character(len=60), intent(in) :: model_type

    !local variables
    character(len=60) :: sub_name
    sub_name = 'allocate_abm_memory'
    call enter_exit(sub_name,1)
    write(*,*) model_type
    call enter_exit(sub_name,2)
end subroutine allocate_abm_memory

end module abm_models
