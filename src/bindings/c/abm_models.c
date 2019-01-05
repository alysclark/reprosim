#include "abm_models.h"

void evaluate_abm_c(const char *model_type);

void evaluate_abm(const char *model_type)
{
    int model_type_len = strlen(model_type);
    evaluate_abm_c(model_type);
}
