#include "read_param.h"

void read_param(void)
{
  ifstream f_in("params/parameters");
  if(!f_in.good())
  {
    cout << "\nCan't open file parameters to read data!\n";
    exit(1);
  }
  
  f_in>>SWEEPS>>THERM>>GAP>>MU>>DT>>READIN>>SEED;
  
  //BETA = (NCOLOR*0.25)/LAMBDA; // Coupling can be absorbed through field redef.
  BETA = 1.0;
  
  TRAJECTORY_LENGTH = (int)(1/DT);
  
  cout << "--------------------------------" << endl ;
  cout << "Number of colors " << NCOLOR <<  "\n";
  cout << "Mass parameter " << MU << "\n";
  cout << "Lattice Coupling " << BETA << "\n";
  cout << "Thermalization sweeps " << THERM << "\n";
  cout << "Number of sweeps " << SWEEPS << "\n";
  cout << "Gap between measurements " << GAP << "\n";
  cout << "Time step in Leapfrog eqs " << DT << "\n";
  cout << "Trajectory length " << TRAJECTORY_LENGTH << "\n";
  cout << "Reading initial config: (1 for yes, 0 for no) " << READIN << "\n";
  
  cout << "--------------------------------" << endl;
  
  setup();
  
  return;
}
