#include <gsl/gsl_rng.h>

void get_gsl_rng_default( const gsl_rng_type* gsl_rng_type_in )
{
  gsl_rng_type_in = gsl_rng_default;
  return;
}
