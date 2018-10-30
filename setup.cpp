#include "setup.h"

void setup(void)
{
  int i,j,k;
  
  // Note our convention: Tr(TaTb) = -delta_ab
  
  cout << "Computing generators for SU(N)\n" << flush;
  
  (void)my_gen();
  
  /*
  // Test orthogonality
  Complex trace;
  for(i=0;i<RANK;i++)
  {
    for(j=0;j<RANK;j++)
    {
      trace = Tr(Lambda[i]*Lambda[j]);
      if(trace.norm()>0.000001)
      {
        cout << "TrT_"<< i << "T_" << j << "= "
              <<  trace << "\n" << flush;
      }
    }
  }*/
  
  if(RANK==NCOLOR*NCOLOR)
  {
    Lambda[RANK-1] = Umatrix();
    
    for(i=0;i<NCOLOR;i++)
    {
      Lambda[RANK-1].set(i,i,(1.0/sqrt(NCOLOR))*Complex(0.0,1.0));
    }
  }
  
  /*
  // structure constants and not-quite-structure constants!
  for(i=0;i<RANK;i++)
    for(j=0;j<RANK;j++)
      for(k=0;k<RANK;k++)
      {
        f[i][j][k] = -Tr(Lambda[i]*(Lambda[j]*Lambda[k]-Lambda[k]*Lambda[j])).real();
        if(fabs(f[i][j][k])>0.0000001)
        {
          cout << "f[" << i << "][" << j << "]["
               << k << "] is " << f[i][j][k] << "\n"
               << flush;
        }
      }
  */
  
  return;
}

