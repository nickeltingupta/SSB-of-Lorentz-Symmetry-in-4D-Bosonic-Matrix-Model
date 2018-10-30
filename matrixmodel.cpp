#include "matrixmodel.h"

int SWEEPS,GAP,THERM,SEED,READIN;
int TRAJECTORY_LENGTH;
double BETA,DT,MASS,MU;
double f[RANK][RANK][RANK];

Umatrix Lambda[RANK];

int main(void)
{
  int i,sweep,sites;
  double WIDTH = 0.1;
  Site_Field X[NSCALAR];
  
  SEED = -1;
  read_param();
  
  if(READIN)
  {
    read_in(X);
  }
  else
  {
    for(i=0;i<NSCALAR;i++)
      X[i] = WIDTH*Site_Field(1);
  }
  
  cout << "Warming up" << "\n" << flush;
  
  DT = 0.5*DT;
  cout << "Using thermalization DT " << DT << "\n" << flush;
  for(sweep=1;sweep<=THERM;sweep++)
  {
    cout << "THERM SWEEP NO. " << sweep << "\n" << flush;
    update(X);
    write_out(X);
  }
  
  cout << "Commencing generation sweeps" << "\n" << flush;
  DT = DT*2;
  cout << "Updated DT is " << DT << "\n" << flush;
  for(sweep=1;sweep<=SWEEPS;sweep++)
  {
    update(X);
    
    // debug: set all matrices in vacuum
#if 0
    Umatrix vac = Umatrix();
    for(int i=0;i<NCOLOR-1;i++)
      vac.set(i,i,Complex(0.0,1.0));
    vac.set(NCOLOR-1,NCOLOR-1,Complex(0.0,1.0 - NCOLOR));
    
    int sites=0,mu;
    Lattice_Vector x;
    
    sites=0;
    while(loop_over_lattice(x,sites))
    {
      for(i=0; i<NSCALAR; i++)
        X[i].set(x,myrandom() * vac);
    }
#endif
    
    // Measure observables
    cout << "GEN SWEEP NO. " << sweep << "\n" << flush;
    if(sweep % GAP == 0)
    {
      measure(X);
    }
  }
  
  write_out(X);
  return 0;
}
