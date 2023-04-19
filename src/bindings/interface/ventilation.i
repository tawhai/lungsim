
%module(package="aether") ventilation
%include symbol_export.h

%{
#include "ventilation.h"
%}

void evaluate_vent(int gdirn, int num_brths, int num_itns, double chest_wall_compliance, double dt, double err_tol, double press_in, double volume_target, double Tinsp, double Texpn, const char *expiration_type);

%include ventilation.h
