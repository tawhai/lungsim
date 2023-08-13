#ifndef AETHER_GAS_EXCHANGE_H
#define AETHER_GAS_EXCHANGE_H

#include "symbol_export.h"

SHO_PUBLIC void steadystate_gasexchange(double cardiac_output, double Vdot_alv, double c_art_o2,
					double c_ven_o2,double p_art_co2,
					double p_art_o2, double p_i_o2, double p_ven_co2,
					double p_ven_o2, double shunt_fraction, double VCO2,
					double VO2);
SHO_PUBLIC void set_perfusion_gradient(int Gdirn, double cardiac_output, double COV, double Qmax, double Qmin);

SHO_PUBLIC void initial_gasexchange(double initial_concentration, double inlet_concentration,
				      double p_ven_o2, double shunt_fraction, double cardiac_output,
				      double surface_area, double V_cap, const char *species);
#endif /* AETHER_GAS_EXCHANGE_H */
