#include "read_in.h"

void read_in(Site_Field X[])
{
  int j,site,LDUM,TDUM;
  double BETADUM,NCOLORDUM,DTDUM;
  Umatrix dummy;
  Lattice_Vector x;
  Site_Field s;
  ifstream f_read;
  
  f_read.open("data/config");
  if(f_read.bad())
    cout << "Error opening config file to read\n" << flush;
  
  f_read >> LDUM >> TDUM >> BETADUM >> DTDUM >> NCOLORDUM;
  
  if((LDUM != L) || (TDUM != T))
    cout << "Wrong size lattice read in - abort\n";
  
  cout << "Config coupling is " << BETADUM << "\n";
  cout << "Time step used for input config " << DTDUM << "\n" << flush;
  cout << "Number of colors " << NCOLORDUM << "\n" << flush;
  
  for(j=0;j<NSCALAR;j++)
  {
    site=0;
    while(loop_over_lattice(x, site))
    {
      f_read >> dummy;
      s.set(x, dummy);
    }
    X[j] = s;
  }
  
  f_read.close();
  return;
}
