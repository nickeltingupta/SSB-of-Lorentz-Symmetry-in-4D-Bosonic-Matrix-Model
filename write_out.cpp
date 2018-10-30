#include "write_out.h"

void write_out(const Site_Field X[])
{
  static ofstream f_write;
  int j,mu,site;
  Lattice_Vector x;
  
  f_write.open("data/dump");
  if (f_write.bad())
    cout << "Error opening config file\n" << flush;
  
  f_write << L << "\t" << T << "\t" << BETA << "\t" << DT << "\t" 
					<< NCOLOR << "\n";
  
 
  for(j=0;j<NSCALAR;j++)
  {
    site=0;
    while(loop_over_lattice(x, site))
      f_write << X[j].get(x) << "\n";
  }
  f_write << "\n" << flush;
  
  f_write.close();
  return;
}
