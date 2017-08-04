/* ----------------------------------------------------------------------------
  INTERP2DPOL
    
    Performs 2-dimensional cubic spline interpolation across a polar grid.
    
    INPUTS
      grid  : input polar grid
      g_rad : radial coordinates of input polar grid
      g_ang : angular coordinates of input polar grid
      i_rad : radii for interpolation
      i_ang : angles for interpolation
      n_rad : number of radial grid points
      n_ang : number of angular grid points
      n_int : number of interpolation points
    
  Mark den Brok
  
  This code is released under a BSD 2-clause license.
  If you use this code for your research, please cite:
  Watkins et al. 2013, MNRAS, 436, 2598
  "Discrete dynamical models of omega Centauri"
  http://adsabs.harvard.edu/abs/2013MNRAS.436.2598W
---------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"
#include "interp.h"


double *interp2dpol( double **grid, double *g_rad, double *g_ang, \
        double *i_rad, double *i_ang, int n_rad, int n_ang, int n_int ) {
    
    int i, j;
    double *interang, *result;
    gsl_spline **splines, *rspline;
    gsl_interp_accel **accs, *acc_r;
    
    
    // allocate memory
    
    accs = (gsl_interp_accel **) malloc( n_rad * sizeof( gsl_interp_accel ) );
    splines = (gsl_spline **) malloc( n_rad * sizeof( gsl_spline ) );
    
    result = (double *) malloc( n_int * sizeof( double ) );
    
    interang = (double *) malloc( n_rad * sizeof( double ) );
    if ( ( result == NULL ) || ( accs == NULL) || ( splines == NULL) \
            || ( interang == NULL ) ) {
        printf( "Malloc failed. Exiting...\n" );
        exit(1);
    }
    
    
    // set up interpolation
    
    acc_r = gsl_interp_accel_alloc();
    rspline = gsl_spline_alloc( gsl_interp_cspline, n_rad );
    
    
    // first do interpolation in the fast part of the array (angles)
    
    for ( i = 0; i < n_rad; i++ ) {
        
        splines[i] = gsl_spline_alloc( gsl_interp_cspline_periodic, n_ang );
        accs[i] = gsl_interp_accel_alloc();
        gsl_spline_init( splines[i], g_ang, grid[i], n_ang );
        
    }
    
    
    // interpolate in radius
    
    for ( j = 0; j < n_int; j++ ) {
        
        // initialize in r-direction by interpolating over all angular values
        for ( i = 0; i < n_rad; i++ ) {
            interang[i] = gsl_spline_eval( splines[i], i_ang[j], accs[i] );
        }
        
        // interpolate over radial values
        gsl_spline_init( rspline, g_rad, interang, n_rad );
        result[j] = gsl_spline_eval( rspline, i_rad[j], acc_r );
        
    }
    
    
    // clean up
    
    for ( i = 0; i < n_rad; i++ ) {
        gsl_spline_free( splines[i] );
        gsl_interp_accel_free( accs[i] );
    }
    
    free( splines );
    free( accs );
    free( interang );
    gsl_spline_free( rspline );
    gsl_interp_accel_free( acc_r );
    
    return result;
    
}
