#include <stdio.h>
#include <math.h>
#include <gsl/gsl_multifit_nlin.h>

#include "rotfit.c"

void temp_rot_fit( int my_n, int my_p, double E_split, double T_init, int low_high_flag, 
		   double my_y[my_n], double my_t[my_n], 
		   double my_c[my_p], double my_cov[my_p][my_p] )
{
  const int low = 0;
  const int high = 1;

  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;

  int status;
  int i, k, iter;
  int n, p;

  n = my_n;
  p = my_p;

  gsl_multifit_function_fdf f;
  gsl_matrix *cov = gsl_matrix_alloc (p, p);

  double y[n], t[n], sigma[n];
  double x_init[p];

  struct data d = { n, p, y, t, sigma };

  x_init[0] = 1.0;
  for (i = 1; i < p; i++)
    {
      x_init[i] = 0.0;
    }
      
  f.f = &rot_f;
  f.df = &rot_df;
  f.fdf = &rot_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  gsl_vector_view x = gsl_vector_view_array (x_init, p);

  for (i = 0; i < n; i++)
    {
      y[i] = my_y[i];
      t[i] = my_t[i];
      sigma[i] = 0.1;
    }

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);

  iter = 0;
  status = GSL_CONTINUE;
  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);

      if(status)
	break;

      status = gsl_multifit_test_delta (s->dx, s->x,
					1e-10, 1e-8);

    }
  while (status == GSL_CONTINUE && iter < 10000);
  /* printf("iterations stopped with status: %s\n", gsl_strerror (status)); */
  /* printf ("iterations reached = %i\n",iter); */

  gsl_multifit_covar (s->J, 0.0, cov);
  
#define FIT(i) gsl_vector_get (s->x, i)
#define COVAR(i,k) gsl_matrix_get (cov, i, k)

  for (i = 0; i <p; i++)
    {
      my_c[i] = FIT(i);
      for (k = 0; k < p; k++)
	{
	  my_cov[k][i] = COVAR(i,k);
	}
    }

  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (cov);

  return;
}

int print_state_rot (int iter, int p, gsl_multifit_fdfsolver * s)
{
  int i;
  
  printf("iter: %i\n\n", iter);

  printf("x:\n");
  for (i = 0; i < p; i++)
    printf("(%i) %f, ", i, gsl_vector_get (s->x, i));
  printf("\n");

  printf("dx:\n");
  for (i = 0; i < p; i++)
    printf("(%i) %f, ", i, gsl_vector_get (s->dx, i));
  printf("\n\n");

  printf("|f(x)| - % 15.8f\n\n\n", gsl_blas_dnrm2 (s->f));
}
