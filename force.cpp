#include "force.h"

void force(const Site_Field X[], Site_Field f_X[])
{
  Lattice_Vector x;
  int sites,i,j,n;
  double td;
  Complex trace;
  Umatrix tmp;
    
  // Scalar force contributions from scalar - scalar self-interactions
  //   -2sum_{i != j} [X_j, [X_i, X_j]]
  for(i=0;i<NSCALAR;i++)
  {
    sites = 0;
    while(loop_over_lattice(x, sites))
    {
      tmp = Umatrix();
      for(j=0;j<NSCALAR;j++)
      {
        if(i==j)
          continue;
        tmp = tmp + comm(X[j].get(x), comm(X[i].get(x), X[j].get(x)));
      }
      
      trace = Tr(tmp);
      if(trace.norm() > TRACETOL)
        tmp = tmp - (1.0 / NCOLOR) * trace * Umatrix(1);
      
      f_X[i].set(x, -2.0 * BETA * tmp);
    }
  }
  
  // Scalar force contributions from mass terms
  td = BETA * MU * MU / 4.5;
  for(i=0;i<NSCALAR;i++)
  {
    sites=0;
    while(loop_over_lattice(x, sites))
    {
      f_X[i].set(x, f_X[i].get(x) - td * X[i].get(x));
    }
  }
  
// Checked: scalar forces are already traceless
  //for(i=0;i<NSCALAR;i++)
  //{
  //  sites=0;
  //  while(loop_over_lattice(x, sites))
  //  {
  //    trace = Tr(f_X[i].get(x));
  //    if (trace.norm() > TRACETOL)
  //      f_X[i].set(x, f_X[i].get(x) - (1.0 / NCOLOR)
  //                   * trace * Umatrix(1));
  //  }
  //}
  
  return;
}
