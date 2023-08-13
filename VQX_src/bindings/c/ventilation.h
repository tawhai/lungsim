#ifndef AETHER_VENTILATION_H
#define AETHER_VENTILATION_H

#include "symbol_export.h"

SHO_PUBLIC void evaluate_vent(int gdirn, int num_brths, int num_itns, double chest_wall_compliance,
			      double dt, double err_tol, double press_in, double volume_target,
			      double Tinsp, double Texpn, const char *expiration_type);
SHO_PUBLIC void evaluate_uniform_flow();
SHO_PUBLIC void two_unit_test();

#endif /* AETHER_VENTILATION_H */
