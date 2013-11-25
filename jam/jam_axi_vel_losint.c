/* -----------------------------------------------------------------------------
  JAM_AXI_VEL_LOSINT
    
    Calculates integrand for line-of-sight integral required for first moment
    calculation.
    
    INPUTS
      zp     : line-of-sight coordinate z' (integration variable)
      params : function parameters passed as a structure
      
    NOTES
      * Based on janis1_jeans_mge_los_integrand IDL code by Michele Cappellari.
    
  Laura L Watkins [lauralwatkins@gmail.com]
  
  This code is released under a BSD 2-clause license.
  If you use this code for your research, please cite:
  Watkins et al. 2013, MNRAS, 436, 2598
  "Discrete dynamical models of omega Centauri"
  http://adsabs.harvard.edu/abs/2013MNRAS.436.2598W
----------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "jam.h"
#include "../mge/mge.h"


double jam_axi_vel_losint( double zp, void *params ) {
    
    struct params_losint *lp;
    struct params_mgeint mp;
    double xp, yp, si, ci, r, z, r2, z2, nu, intg, result, error, *knu;
    int i;
    
    // get parameters
    lp = params;
    xp = lp->xp;
    yp = lp->yp;
    
    // intrinsic R and z
    si = sin( lp->incl );
    ci = cos( lp->incl );
    r = sqrt( pow( zp * si - yp * ci, 2 ) + pow( xp, 2 ) );         // eqn 25
    z = sqrt( pow( zp * ci + yp * si, 2 ) );
    
    // do some prep for the integrand to avoid repeat calculations
    r2 = r * r;
    z2 = z * z;
    
    // kappa^2 * nu
    knu = (double *) malloc( lp->lum->ntotal * sizeof( double ) );
    for ( i = 0; i < lp->lum->ntotal; i++ ) {
        if ( lp->kappa[i] == 0. ) knu[i] = 0.;
        else knu[i] = pow( lp->kappa[i], 3 ) / fabs( lp->kappa[i] ) \
            * lp->lum->area[i] \
            * exp( -0.5 / lp->s2l[i] * ( r2 + z2 / lp->q2l[i] ) );
    }
    
    // parameters for integrand function
    mp.r2 = r2;
    mp.z2 = z2;
    mp.lum = lp->lum;
    mp.pot = lp->pot;
    mp.bani = lp->bani;
    mp.s2l = lp->s2l;
    mp.q2l = lp->q2l;
    mp.s2q2l = lp->s2q2l;
    mp.s2p = lp->s2p;
    mp.e2p = lp->e2p;
    mp.knu = knu;
    
    // perform integration
    gsl_integration_workspace *w = gsl_integration_workspace_alloc( 1000 );
    gsl_function F;
    F.function = &jam_axi_vel_mgeint;
    F.params = &mp;
    gsl_integration_qag( &F, 0., 1., 0., 1e-5, 1000, 2, w, &result, &error );
    gsl_integration_workspace_free( w );
    
    // mge volume density
    nu = mge_dens( lp->lum, r, z );
    
    // keep track of kappa signs - see note 8 p77 of Cappellari 2008
    intg = nu * result / fabs( nu * result ) * sqrt( fabs( nu * result ) );
    
    free( knu );
    
    intg *= pow( zp, lp->zpow );
    
    return intg;
    
}
