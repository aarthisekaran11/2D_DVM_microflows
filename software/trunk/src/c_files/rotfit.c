#include <stdio.h>
#include <math.h>
#include <gsl/gsl_multifit_nlin.h>

struct data {
  int n;
  int p;
  double * y;
  double * t;
  double * sigma;
};

int rot_f( const gsl_vector *x, void *data, gsl_vector * f )
{
  int n = ((struct data *) data)->n;
  int p = ((struct data *) data)->p;
  double *y = ((struct data *) data)->y;
  double *t = ((struct data *) data)->t;
  double *sigma = ((struct data *) data)->sigma;

  int i, k;

  double param[p];
  for (i = 0; i < p; i++)
    {
      param[i] = gsl_vector_get (x, i);
    }

  for (i = 0; i < n; i++)
    {
      double s = 1.0/sigma[i];
      double Yi = 0.0;

      for (k = 0; k < p; k++)
	{
	  double power = (double) k;
	  double p_term = pow( t[i], power );
	  Yi = Yi + param[k] * p_term;
	}

      gsl_vector_set (f, i, (Yi - y[i]) * s);

    }

  return GSL_SUCCESS;
}

int rot_df( const gsl_vector * x, void *data, gsl_matrix * J )
{
  int n = ((struct data *) data)->n;
  int p = ((struct data *) data)->p;
  double *y = ((struct data *) data)->y;
  double *t = ((struct data *) data)->t;
  double *sigma = ((struct data *) data)->sigma;

  int i, k;

  double param[p];
  for (i = 0; i < p; i++)
    {
      param[i] = gsl_vector_get (x, i);
    }
  
  /* Loop over all data points */
  for (i = 0; i < n; i++)
    {
      double s = 1.0/sigma[i];
      
      /* The Jacobian is J(i,j) = dfi / dxj
	 where fi = (Yi - yi)/sigma[i] and 
	 where Yi is defined as above and 
	 where xj are the parameters for the fit */
      //      gsl_matrix_set (J, i, 0, s);

      if (p > 0)
	{
	  for (k = 0; k < p; k++)
	    {
	      double power = (double) k;
	      double p_term = pow( t[i], power );
	      gsl_matrix_set (J, i, k, p_term * s);
	    }
	}

    }

  return GSL_SUCCESS;
}

int rot_fdf( const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J )
{
  rot_f (x, data, f);
  rot_df (x, data, J);

  return GSL_SUCCESS;
}

