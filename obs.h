#include "utilities.h"
#include "action.h"

void obs(const Site_Field X[NSCALAR], double &act_s);

extern "C" void zgeev_( char*, char*, int*, double at[], int *,
                       double b[], double dummy[], int *,
                       double dummy2[], int*, double work[],
                       int *, double work2[], int *);
