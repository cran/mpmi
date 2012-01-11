#include <R.h>

// Value of missing integers (INT_MIN)
int F77_SUB(rnaint)(void) { return(R_NaInt); }

// Function to check for missing real values
int F77_SUB(rfinite)(double *x) { return(R_FINITE(*x)); }

