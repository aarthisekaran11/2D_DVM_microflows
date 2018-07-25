#include <stdio.h>
#include <math.h>
#include <gsl/gsl_multifit_nlin.h>

struct ln_data {
  int n;
  int p;
  double c_ln;
  double * y;
  double * t;
  double * sigma;
};

int ln_f( const gsl_vector * x, void *ln_data, gsl_vector * f )
{
  int n = ((struct ln_data *) ln_data)->n;
  int p = ((struct ln_data *) ln_data)->p;
  double c_ln = ((struct ln_data *) ln_data)->c_ln;
  double *y = ((struct ln_data *) ln_data)->y;
  double *t = ((struct ln_data *) ln_data)->t;
  double *sigma = ((struct ln_data *) ln_data)->sigma;

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
	 y = A*exp(-Bt) + C + C1*t + C2*t^2 + ... */
      double ln = log(c_ln/t[i]+1.0);
      double s = sigma[i];
      double Yi = param[0]/ln;
      if (p > 1)
	{
	  for (k = 1; k < p; k++)
	    {
	      double power = (double) k;
	      double p_term = pow( t[i], power );
	      Yi = Yi + param[k] * p_term;
	    }
	}
      gsl_vector_set (f, i, (Yi - y[i])/s);

    }
  return GSL_SUCCESS;
}

int ln_df( const gsl_vector * x, void *ln_data, gsl_matrix * J )
{
  int n = ((struct ln_data *) ln_data)->n;
  int p = ((struct ln_data *) ln_data)->p;
  double c_ln = ((struct ln_data *) ln_data)->c_ln;
  double *y = ((struct ln_data *) ln_data)->y;
  double *t = ((struct ln_data *) ln_data)->t;
  double *sigma = ((struct ln_data *) ln_data)->sigma;

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
      double ln = log(c_ln/t[i]+1.0);
      double s = sigma[i];
      
      /* The Jacobian is J(i,j) = dfi / dxj
	 where fi = (Yi - yi)/sigma[i] and 
	 where Yi is defined as above and 
	 where xj are the parameters for the fit */
      gsl_matrix_set (J, i, 0, 1 / (ln * s) );

      if (p > 1)
	{
	  for (k = 1; k < p; k++)
	    {
	      double power = (double) k;
	      double p_term = pow( t[i], power );
	      gsl_matrix_set (J, i, k, p_term / s);
	    }
	}

    }
  return GSL_SUCCESS;
}

int ln_fdf( const gsl_vector * x, void *ln_data, gsl_vector * f, gsl_matrix * J )
{
  ln_f (x, ln_data, f);
  ln_df (x, ln_data, J);

  return GSL_SUCCESS;
}
