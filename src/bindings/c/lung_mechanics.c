
#include "lung_mechanics.h"
#include "string.h"

#include "stdio.h"

void deform_tissue_in_cavity_c(int *nsteps, const char *posture, int *posture_len, const char *filename, int *filename_len);


void deform_tissue_in_cavity(int nsteps, const char *posture, const char *filename)
{
  int posture_len = strlen(posture);
  int filename_len = strlen(filename);
  deform_tissue_in_cavity_c( &nsteps, posture, &posture_len, filename, &filename_len);
}

