#include "evolve_fields.h"

void evolve_fields(Site_Field X[], Site_Field p_X[], Site_Field f_X[])
{
  Site_Field new_f_X[NSCALAR];
  int i;
  
	// Update fields 
  for (i=0;i<NSCALAR;i++)
    X[i] = X[i] + DT * p_X[i] + 0.5 * DT * DT * f_X[i];
  
  // Update forces
  force(X, new_f_X);
  
  // Update momenta 
  for (i=0;i<NSCALAR;i++)
    p_X[i] = p_X[i] + 0.5 * DT * (new_f_X[i] + f_X[i]);
  
  // Store final forces for next iteration
  for (i=0;i<NSCALAR;i++)
    f_X[i] = new_f_X[i];

  return;
}
