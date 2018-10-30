#include "action.h"

double action(const Site_Field X[NSCALAR])
{
  int i,j,site;
  double act_s=0.0, tr_double=0.0;
  double m_coeff=0.0, mass=0.0;
  Lattice_Vector x;
  Umatrix COMM;
  
  // Scalar - scalar potential -sum_{i<j} Tr[X_i, X_j]^2
  for(i=0;i<NSCALAR;i++)
  {
    for(j=i+1;j<NSCALAR;j++)
    {
      site = 0;
      while(loop_over_lattice(x, site))
      {
        COMM = comm(X[i].get(x), X[j].get(x));
        act_s = act_s - BETA * Tr(COMM * COMM).real();
      }
    }
  }
  
  // Mass terms
  m_coeff = BETA * MU * MU / 9.0;
  for(i=0;i<NSCALAR;i++)
  {
    site = 0;
    while(loop_over_lattice(x, site))
    {
      tr_double = m_coeff * Tr(X[i].get(x) * X[i].get(x)).real();
      act_s = act_s - tr_double;
      mass = mass - tr_double;
    }
  }
  
  return (act_s);
}
