#include "update.h"

void update(Site_Field X[NSCALAR])
{
  static int first_time=1, accept=0, no_calls=0;
  int i;
  static double S_old=0.0, acc=0.0;
  double H_old, K_old, K_new, S_new, H_new, hmc_test;
  static ofstream f_hmc,f_acc;
  
  static Site_Field f_X[NSCALAR];
  Site_Field p_X[NSCALAR], old_X[NSCALAR], old_f_X[NSCALAR];
  
  no_calls++;
  
  // Refresh momenta 
  for(i=0;i<NSCALAR;i++)
    p_X[i] = Site_Field(1);
  
  K_old = kinetic_energy(p_X);
  
  if(first_time)
  {
    f_hmc.open("data/hmc_test");
    if(f_hmc.bad())
      cout << "Failed to open hmc_test file\n" << flush;
    f_acc.open("data/acc_rate");
    if(f_acc.bad())
      cout << "Failed to open acceptance rate file\n" << flush;
    
    S_old = action(X);
    //cout << "Computed action\n" << flush;
    
    force(X,f_X);
    first_time = 0;
  }
  //cout << "S_old is " << S_old << "\n";
  
  if((no_calls % 25 == 0) && (!first_time))
  {
    acc = (double)accept / (double)no_calls;
    
    cout << "Acceptance rate " << acc << endl;
    
    f_acc << acc << endl;
    
    no_calls = 0;
    accept = 0;
  }
  H_old = S_old + K_old;
  
  // Save starting fields
  for(i=0;i<NSCALAR;i++)
    old_X[i] = X[i];
  
  for(i=0;i<NSCALAR;i++)
    old_f_X[i] = f_X[i];
  
  cout << "Action at beginning traj " << S_old << "\n" << flush;
  cout << "H at beginning traj " << H_old << "\n" << flush;
  
  // MD evolution
  for(i=0;i<TRAJECTORY_LENGTH;i++)
    evolve_fields(X,p_X,f_X);
  
  S_new = action(X);
  K_new = kinetic_energy(p_X);
  H_new = S_new + K_new;
  
  cout << "Action at end of traj " << S_new << "\n" << flush;
  cout << "H at end of traj " << H_new << "\n" << flush;
  
  // Metropolis test
  hmc_test = exp(H_old - H_new);
  f_hmc << hmc_test << "\n" << flush;
  
  if(myrandom() < hmc_test)
  {
    S_old = S_new;
    accept++;
    cout << "ACCEPTED with deltaH of " << H_new - H_old << endl;
    return;
  }
  else
  {    // Restore starting fields
    cout << "hmc_test " << hmc_test << " failed\n" << flush;
    cout << "REJECTED with deltaH of " << H_new - H_old << endl;
    for(i=0;i<NSCALAR;i++)
    {
      X[i] = old_X[i];
      f_X[i] = old_f_X[i];
    }
    return;
  }
  return;
}
