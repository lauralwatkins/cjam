/* ----------------------------------------------------------------------------
  JAM_AXI_RMS_WMMT
    
    Calculates weighted second moment.
    
    INPUTS
      xp    : projected x' [pc]
      yp    : projected y' [pc]
      nxy   : number of x' and y' values given
      incl  : inclination [radians]
      lum   : projected luminous MGE
      pot   : projected potential MGE
      beta  : velocity anisotropy (1 - vz^2 / vr^2)
      vv    : velocity integral selector (1=xx, 2=yy, 3=zz, 4=xy, 5=xz, 6=yz)
    
    NOTES
    * Based on janis2_weighted_second_moment_squared IDL code by Michele
      Cappellari.
    
  Mark den Brok
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


double *jam_axi_rms_wmmt( double *xp, double *yp, int nxy, double incl, \
        struct multigaussexp *lum, struct multigaussexp *pot, \
        double *beta, int vv ) {
    
    struct params_rmsint p;
    struct multigaussexp ilum, ipot;
    double ci, si, *kani, *s2l, *q2l, *s2q2l, *s2p, *e2p;
    double result, error, *sb_mu2;
    int i;
    
    // convert from projected MGEs to intrinsic MGEs
    ilum = mge_deproject( lum, incl );
    ipot = mge_deproject( pot, incl );
    
    // angles
    ci = cos( incl );
    si = sin( incl );
    
    
    // mge component combinations
    
    kani = (double *) malloc ( lum->ntotal * sizeof( double ) );
    s2l = (double *) malloc( lum->ntotal * sizeof( double ) );
    q2l = (double *) malloc( lum->ntotal * sizeof( double ) );
    s2q2l = (double *) malloc( lum->ntotal * sizeof( double ) );
    s2p = (double *) malloc( pot->ntotal * sizeof( double ) );
    e2p = (double *) malloc( pot->ntotal * sizeof( double ) );
    
    for ( i = 0; i < ilum.ntotal; i++ ) {
        kani[i] = 1. / ( 1. - beta[i] );
        s2l[i] = pow( ilum.sigma[i], 2. );
        q2l[i] = pow( ilum.q[i], 2. );
        s2q2l[i] = s2l[i] * q2l[i];
    }
    
    for ( i = 0; i < ipot.ntotal; i++ ) {
        s2p[i] = pow( ipot.sigma[i], 2. );
        e2p[i] = 1. - pow( ipot.q[i], 2. );
    }
    
    // parameters for the integrand function
    p.ci2 = ci * ci;
    p.si2 = si * si;
    p.cisi = ci * si;
    p.lum = &ilum;
    p.pot = &ipot;
    p.kani = kani;
    p.s2l = s2l;
    p.q2l = q2l;
    p.s2q2l = s2q2l;
    p.s2p = s2p;
    p.e2p = e2p;
    p.vv = vv;
    
    
    // perform integration
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc( 1000 );
    gsl_function F;
    F.function = &jam_axi_rms_mgeint;
    
    sb_mu2 = (double *) malloc( nxy * sizeof( double ) );
    for ( i = 0; i < nxy; i++ ) {
        p.x2 = xp[i] * xp[i];
        p.y2 = yp[i] * yp[i];
        p.xy = xp[i] * yp[i];
        F.params = &p;
        gsl_integration_qag( &F, 0., 1., 0., 1e-5, 1000, 2, w, \
            &result, &error );
        sb_mu2[i] = result;
    }
    
    gsl_integration_workspace_free( w );
    
    free( kani );
    free( s2l );
    free( q2l );
    free( s2q2l );
    free( s2p );
    free( e2p );
    
    return sb_mu2;
    
}
