#include "obs.h"

void obs(const Site_Field X[],
         double &act_s)
{
  Lattice_Vector x;
  int site,i,j;
  
	act_s = 0.0;
  
  // Scalar action
  for(i=0;i<NSCALAR;i++)
  {
    for(j=i+1;j<NSCALAR;j++)
    {
      site=0;
      while(loop_over_lattice(x,site))
      {
        act_s = act_s - BETA*Tr(comm(X[i].get(x),X[j].get(x))*
                            comm(X[i].get(x),X[j].get(x))).real();
      }
    }
  }
  
  // Mass terms
  for(i=0;i<NSCALAR;i++)
  {
    site=0;
    while(loop_over_lattice(x,site))
    {
      act_s = act_s - BETA*MU*MU*Tr(X[i].get(x)*X[i].get(x)).real();
    }
  }
  
  return;
}
