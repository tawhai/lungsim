
#include "gas_exchange.h"
#include "string.h"

void steadystate_gasexchange_c(double *cardiac_output, double *Vdot_alv, double *c_art_o2,
			       double *c_ven_o2,double *p_art_co2, 
			       double *p_art_o2, double *p_i_o2, double *p_ven_co2,
			       double *p_ven_o2, double *shunt_fraction,
			       double *VCO2, double *VO2);
void set_perfusion_gradient_c(int *Gdirn, double *cardiac_output, double *COV, double *Qmax, double *Qmin);
void initial_gasexchange_c(double *initial_concentration, double *inlet_concentration,
			   double *p_ven_o2, double *shunt_fraction, double *cardiac_output,
			   double *surface_area, double *V_cap, const char *species, int *species_len);

  
void steadystate_gasexchange(double cardiac_output, double Vdot_alv, double c_art_o2,
			     double c_ven_o2,double p_art_co2, double p_art_o2, 
			     double p_i_o2, double p_ven_co2, double p_ven_o2, double shunt_fraction,
			     double VCO2, double VO2)
{
  steadystate_gasexchange_c(&cardiac_output, &Vdot_alv, &c_art_o2, &c_ven_o2, &p_art_co2, &p_art_o2, &p_i_o2, &p_ven_co2, &p_ven_o2, &shunt_fraction, &VCO2, &VO2);
}

void set_perfusion_gradient(int Gdirn, double cardiac_output, double COV, double Qmax, double Qmin)
{
  set_perfusion_gradient_c(&Gdirn, &cardiac_output, &COV, &Qmax, &Qmin);
}

void initial_gasexchange(double initial_concentration, double inlet_concentration,
			   double p_ven_o2, double shunt_fraction, double cardiac_output,
			   double surface_area, double V_cap, const char *species)
{
  int species_len = strlen(species);
  initial_gasexchange_c(&initial_concentration, &inlet_concentration, &p_ven_o2, &shunt_fraction, &cardiac_output, &surface_area, &V_cap, species, &species_len);
}
