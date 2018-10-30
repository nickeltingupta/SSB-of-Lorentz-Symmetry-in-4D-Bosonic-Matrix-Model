#include "measure.h"

void measure(const Site_Field X[])
{
  static int first_time = 1;
  static ofstream f_act, f_ssq, f_J;
  double ssq_val=0.0, scalar_dum, act_s;
  double Iij=0.0, Iii=0.0, J=0.0;
  int site,i,j,k;
  Lattice_Vector x;
 
  if(first_time)
  {
    f_act.open("data/action");
    if(f_act.bad())
    {
      cout << "Failed to open action file\n" << flush;
    }

    f_ssq.open("data/ssq");
    if(f_ssq.bad())
    {
      cout << "Failed to open scalar squares file\n" << flush;
    }
    f_J.open("data/jvals");
    if(f_J.bad())
    {
      cout << "Failed to open J obsevable file\n" << flush;
    }
    first_time=0;
  }
  
  obs(X,act_s);

  f_act << act_s << endl;

  for(i=0;i<NSCALAR;i++)
  {
    site=0;
    while(loop_over_lattice(x, site))
    {
      scalar_dum = Tr(X[i].get(x) * X[i].get(x)).real();
      ssq_val = ssq_val - scalar_dum;
    }
		f_ssq << ssq_val << "\t";
		ssq_val = 0.0;
  }
  f_ssq << endl;
  
  for(i=0;i<NSCALAR;i++)
  {
    for(j=0;j<NSCALAR;j++)
    {
      site=0;
      while(loop_over_lattice(x, site))
      {
        Iij = (1.0/NCOLOR)*(Tr(X[i].get(x) * X[j].get(x)).real());
        Iii = (1.0/NCOLOR)*(Tr(X[i].get(x) * X[i].get(x)).real());
        J = J + (0.25*Iij*Iij - (0.25*Iii)*(0.25*Iii));
      }
    }
  }
  f_J << 1.0*J/NCOLOR << endl;
  
  return;
}
