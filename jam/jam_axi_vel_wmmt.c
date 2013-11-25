/* ----------------------------------------------------------------------------
  JAM_AXI_VEL_WMMT
    
    Calculates weighted first moments.
    
    INPUTS
      xp    : projected x' [pc]
      yp    : projected y' [pc]
      nxy   : number of x' and y' values given
      incl  : inclination [radians]
      lum   : projected luminous MGE
      pot   : projected potential MGE
      beta  : velocity anisotropy (1 - vz^2 / vr^2)
      kappa : rotation parameter
    
    NOTES
      * Based on janis1_weighted_first_moment IDL code by Michele Cappellari.
    
  Laura L Watkins [lauralwatkins@gmail.com]
  
  This code is released under a BSD 2-clause license.
  If you use this code for your research, please cite:
  Watkins et al. 2013, MNRAS, 436, 2598
  "Discrete dynamical models of omega Centauri"
  http://adsabs.harvard.edu/abs/2013MNRAS.436.2598W
---------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "jam.h"
#include "../mge/mge.h"
#include "../tools/tools.h"


double** jam_axi_vel_wmmt( double *xp, double *yp, int nxy, double incl, \
        struct multigaussexp *lum, struct multigaussexp *pot, \
        double *beta, double *kappa ) {
    
    struct params_losint lp;
    struct multigaussexp ilum, ipot;
    double *bani, *s2l, *q2l, *s2q2l, *s2p, *e2p;
    double *iz0, *iz1;
    double lim, result, error, si, ci, trpig, **sb_mu1;
    int i;
    
    // ---------------------------------
    
    // MGE parameters for inner integrand function
    
    // convert from projected MGEs to intrinsic MGEs
    ilum = mge_deproject( lum, incl );
    ipot = mge_deproject( pot, incl );
    
    // calculate useful quantities
    bani = (double *) malloc( ilum.ntotal * sizeof( double ) );
    s2l = (double *) malloc( ilum.ntotal * sizeof( double ) );
    q2l = (double *) malloc( ilum.ntotal * sizeof( double ) );
    s2q2l = (double *) malloc( ilum.ntotal * sizeof( double ) );
    s2p = (double *) malloc( ipot.ntotal * sizeof( double ) );
    e2p = (double *) malloc( ipot.ntotal * sizeof( double ) );
    
    for ( i = 0; i < ilum.ntotal; i++ ) {
        bani[i] = 1. / ( 1. - beta[i] );
        s2l[i] = pow( ilum.sigma[i], 2 );
        q2l[i] = pow( ilum.q[i], 2 );
        s2q2l[i] = s2l[i] * q2l[i];
    }
    
    for ( i = 0; i < ipot.ntotal; i++ ) {
        s2p[i] = pow( ipot.sigma[i], 2 );
        e2p[i] = 1. - pow( ipot.q[i], 2 );
    }
    
    // parameters for integrand function
    lp.incl = incl;
    lp.lum = &ilum;
    lp.pot = &ipot;
    lp.bani = bani;
    lp.s2l = s2l;
    lp.q2l = q2l;
    lp.s2q2l = s2q2l;
    lp.s2p = s2p;
    lp.e2p = e2p;
    lp.kappa = kappa;
    
    // ---------------------------------
    
    // set up integration
    gsl_integration_workspace *w = gsl_integration_workspace_alloc( 1000 );
    gsl_function F;
    F.function = &jam_axi_vel_losint;
    
    // trig angles
    si = sin( incl );
    ci = cos( incl );
    
    // outer limit of integration
    lim = 4. * maximum( ilum.sigma, ilum.ntotal );
    
    iz0 = (double *) malloc( nxy * sizeof( double ) );
    iz1 = (double *) malloc( nxy * sizeof( double ) );
    for ( i = 0; i < nxy; i++ ) {
        
        // parameters for integrand function
        lp.xp = xp[i];
        lp.yp = yp[i];
        
        // do z^0 integral
        lp.zpow = 0.;
        F.function = &jam_axi_vel_losint;
        F.params = &lp;
        gsl_integration_qag( &F, -lim, lim, 1e-3, 1e-5, 1000, 2, w, \
            &result, &error );
        iz0[i] = result;
        
        // do z^1 integral
        lp.zpow = 1.;
        F.function = &jam_axi_vel_losint;
        F.params = &lp;
        gsl_integration_qag( &F, -lim, lim, 1e0, 1., 1000, 2, w, \
            &result, &error );
        if ( abs( result ) != 0. ) {
            lp.zpow = 1.;
            F.function = &jam_axi_vel_losint;
            F.params = &lp;
            gsl_integration_qag( &F, -lim, lim, 1e-3, 1e-5, 1000, 2, w, \
                &result, &error );
        }
        iz1[i] = result;
        
    }
    
    // tidy up
    gsl_integration_workspace_free( w );
    
    // ---------------------------------
    
    sb_mu1 = (double **) malloc( nxy * sizeof( double* ) );
    for ( i = 0; i < nxy; i++ ) \
        sb_mu1[i] = (double *) malloc( 3 * sizeof( double ) );
    
    // calculate for each velocity component
    trpig = 2. * sqrt( M_PI * G );
    for ( i = 0; i < nxy; i++ ) {
        sb_mu1[i][0] = trpig * ( yp[i] * ci * iz0[i] - si * iz1[i] );
        sb_mu1[i][1] = -trpig * xp[i] * ci * iz0[i];
        sb_mu1[i][2] = trpig * xp[i] * si * iz0[i];
    }
    
    // ---------------------------------
    
    free( bani );
    free( s2l );
    free( q2l );
    free( s2q2l );
    free( s2p );
    free( e2p );
    
    free( iz0 );
    free( iz1 );
    
    return sb_mu1;
    
}
