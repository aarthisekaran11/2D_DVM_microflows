#include <stdio.h>
#include <gsl/gsl_multifit.h>

void curve_fit( int my_n, int my_p, double my_X[my_p][my_n], double my_y[my_n], double my_c[my_p], 
		double my_cov[my_p][my_p], double my_chisq )
{
  int i, j;
  int n, p;
  double chisq;
  gsl_matrix *X, *cov;
  gsl_vector *y, *c;

  n = my_n;
  p = my_p;

  X = gsl_matrix_alloc (n, p);
  y = gsl_vector_alloc (n);

  c   = gsl_vector_alloc (p);
  cov = gsl_matrix_alloc (p, p);

  for (i = 0; i < n; i++)
    {
      for (j = 0; j < p; j++)
  	{
  	  gsl_matrix_set (X, i, j, my_X[j][i]);
  	}
      gsl_vector_set (y, i, my_y[i+1]);
    }
  
  gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, p);
  gsl_multifit_linear (X, y, c, cov, &chisq, work);
  gsl_multifit_linear_free (work);

#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

  for (i = 0; i < p; i++)
    {
      for (j = 0; j < p; j++)
  	{
	  my_cov[j][i] = COV(i,j);
  	}
      my_c[i] = C(i);
    }

  my_chisq = chisq;

  gsl_matrix_free (X);
  gsl_vector_free (y);

  return;
}
