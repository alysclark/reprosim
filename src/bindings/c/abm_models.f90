module abm_models_c
implicit none
private

contains

!!!###################################################################################

subroutine evaluate_abm_c(model_type,model_type_len) bind(C, name="evaluate_abm_c")

use iso_c_binding, only: c_ptr
use utils_c, only: strncpy
use other_consts, only: MAX_STRING_LEN
use arrays, only: dp
use abm_models, only: evaluate_abm
implicit none

type(c_ptr), value, intent(in) :: model_type
integer,intent(in) :: model_type_len
character(len=MAX_STRING_LEN) :: model_type_f

call strncpy(model_type_f, model_type, model_type_len)

#if defined _WIN32 && defined __INTEL_COMPILER
call so_evaluate_abm(model_type_f)
#else
call evaluate_abm(model_type_f)
#endif

end subroutine evaluate_abm_c

!###################################################################################
end module abm_models_c
