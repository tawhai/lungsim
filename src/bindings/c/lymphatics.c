#include "lymphatics.h"
#include "string.h"

void alveolar_flux_c(double *dt, double *time, double *T_interval);
void lymphatic_transport_c(const char *filename, int *filename_len);


void alveolar_flux(double dt, double time, double T_interval)
{
  alveolar_flux_c(&dt, &time, &T_interval);
}

void lymphatic_transport(const char *filename)
{
  int filename_len = (int)strlen(filename);
  lymphatic_transport_c(filename, &filename_len);
}
