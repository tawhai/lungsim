
#include "string.h"
#include "ventilation.h"

void evaluate_vent_c(int *gdirn, int *num_brths, int *num_itns, double *chest_wall_compliance, double *dt, double *err_tol, double *press_in, double *volume_target, double *Tinsp, double *Texpn, const char *expiration_type, int *expiration_type_len);
void evaluate_uniform_flow_c();
void two_unit_test_c();

void evaluate_vent(int gdirn, int num_brths, int num_itns, double chest_wall_compliance, double dt, double err_tol, double press_in, double volume_target, double Tinsp, double Texpn, const char *expiration_type)
{
  int expiration_type_len = (int)strlen(expiration_type);
  evaluate_vent_c(&gdirn, &num_brths, &num_itns, &chest_wall_compliance, &dt, &err_tol, &press_in, &volume_target, &Tinsp, &Texpn, expiration_type, &expiration_type_len);
}

void evaluate_uniform_flow()
{
  evaluate_uniform_flow_c();
}

void two_unit_test()
{
  two_unit_test_c();
}
