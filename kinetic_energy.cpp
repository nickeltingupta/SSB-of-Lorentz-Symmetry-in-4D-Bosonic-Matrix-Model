#include "kinetic_energy.h"

double kinetic_energy(const Site_Field p_X[])
{
  int i;
  Complex dum = Complex();
  
  for(i=0;i<NSCALAR;i++)
    dum = dum + 0.5 * Tr(Adj(p_X[i]) * p_X[i]);
  
  return dum.real();
}
