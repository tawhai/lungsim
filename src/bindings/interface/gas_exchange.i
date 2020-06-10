
%module(package="aether") gas_exchange
%include symbol_export.h
%include gas_exchange.h

%{
#include "gas_exchange.h"
%}

void initial_gasexchange(double initial_concentration, double inlet_concentration,
			 double p_ven_o2, double shunt_fraction, double cardiac_output,
			 double surface_area, double V_cap, const char *species);

%include gas_exchange.h
