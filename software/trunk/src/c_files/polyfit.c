#include <stdio.h>
#include <math.h>
#include <gsl/gsl_multifit_nlin.h>

struct poly_data {
  int n;
  int p;
  double * y;
  double * t;
  double * sigma;
};

int poly_f( const gsl_vector * x, void *poly_data, gsl_vector * f )
{
  int n = ((struct poly_data *) poly_data)->n;
  int p = ((struct poly_data *) poly_data)->p;
  double *y = ((struct poly_data *) poly_data)->y;
  double *t = ((struct poly_data *) poly_data)->t;
  double *sigma = ((struct poly_data *) poly_data)->sigma;

  int i, k;

  /* Read parameters from x gsl_vector */
  double param[p];
  for (i = 0; i < p; i++)
    {
      param[i] = gsl_vector_get (x, i);
    }
  
  /* Loop over all data points */
  for (i = 0; i < n; i++)
    {
      /* The current curve fit is of the form:
	 y = C0 + C1*t + C2*t^2 + ... */
      double s = sigma[i];
      double Yi = 0.0;
      for (k = 0; k < p; k++)
	{
	  double power = (double) k;
	  double p_term = pow( t[i], power );
	  Yi = Yi + param[k] * p_term;
	}
      gsl_vector_set (f, i, (Yi - y[i])/s);

    }
  return GSL_SUCCESS;
}

int poly_df( const gsl_vector * x, void *poly_data, gsl_matrix * J )
{
  int n = ((struct poly_data *) poly_data)->n;
  int p = ((struct poly_data *) poly_data)->p;
  double *y = ((struct poly_data *) poly_data)->y;
  double *t = ((struct poly_data *) poly_data)->t;
  double *sigma = ((struct poly_data *) poly_data)->sigma;

  int i, k;

  /* Read parameters from x gsl_vector */
  double param[p];
  for (i = 0; i < p; i++)
    {
      param[i] = gsl_vector_get (x, i);
    }
  
  /* Loop over all data points */
  for (i = 0; i < n; i++)
    {
      double s = sigma[i];
      
      /* The Jacobian is J(i,j) = dfi / dxj
	 where fi = (Yi - yi)/sigma[i] and 
	 where Yi is defined as above and 
	 where xj are the parameters for the fit */
      for (k = 0; k < p; k++)
	{
	  double power = (double) k;
	  double p_term = pow( t[i], power );
	  gsl_matrix_set (J, i, k, p_term / s);
	}

    }
  return GSL_SUCCESS;
}

int poly_fdf( const gsl_vector * x, void *poly_data, gsl_vector * f, gsl_matrix * J )
{
  poly_f (x, poly_data, f);
  poly_df (x, poly_data, J);

  return GSL_SUCCESS;
}
