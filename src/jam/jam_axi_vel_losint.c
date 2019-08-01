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
#include <gsl/gsl_errno.h>
#include "jam.h"
#include "../mge/mge.h"


double jam_axi_vel_losint(double zp, void *params) {
    
    struct params_losint *lp;
    struct params_mgeint mp;
    double xp, yp, si, ci, r, z, r2, z2, nu, intg, result, error, nu_i;
    double sign_kappa, sum;
    int i;
    int status = 0;
    
    // get parameters
    lp = params;
    xp = lp->xp;
    yp = lp->yp;
    
    // intrinsic R and z
    si = sin(lp->incl);
    ci = cos(lp->incl);
    r = sqrt(pow(zp*si - yp*ci, 2) + pow(xp, 2));                   // eqn 25
    z = sqrt(pow(zp*ci + yp*si, 2));
    
    // do some prep for the integrand to avoid repeat calculations
    r2 = r * r;
    z2 = z * z;
    
    // parameters for integrand function
    mp.r2 = r2;
    mp.z2 = z2;
    mp.pot = lp->pot;
    mp.s2p = lp->s2p;
    mp.e2p = lp->e2p;
    
    // perform integration
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function = &jam_axi_vel_mgeint;
    gsl_set_error_handler_off(); 
    sum = 0.;
    for (i=0; i<lp->lum->ntotal; i++) {
        if (lp->kappa[i]==0.) sign_kappa = 0.;
        else sign_kappa = lp->kappa[i]/fabs(lp->kappa[i]);
        nu_i = lp->lum->area[i] * exp(-0.5/lp->s2l[i]*(r2+z2/lp->q2l[i]));
        
        mp.bani = lp->bani[i];
        mp.s2l = lp->s2l[i];
        mp.q2l = lp->q2l[i];
        mp.s2q2l = lp->s2q2l[i];
        F.params = &mp;
        status += gsl_integration_qag(&F, 0., 1., 0., 1e-5, 1000, 6, w, &result, &error);
        sum += sign_kappa * pow(lp->kappa[i], 2) * nu_i * fabs(result);
    }
    
    gsl_integration_workspace_free(w);
    
    // mge volume density
    nu = mge_dens(lp->lum, r, z);
    
    // keep track of kappa signs - see note 8 p77 of Cappellari 2008
    intg = nu*sum/fabs(nu*sum)*sqrt(fabs(nu*sum));
    
    intg *= pow(zp, lp->zpow);
    
    /*// if you want to test the gsl flagging thing uncomment these two lines
    int test = 1;// Tadeja
    status += test; // Tadeja
    */
    if (status > 0 && *lp->gslFlag_losint<=10) { // the limiting value is very arbitrary so that the values do not exceed allowed size for integers
        *lp->gslFlag_losint += 1;   }

    if (status > 0 && *lp->gslFlag_losint>10) { // the limiting value is very arbitrary so that the values do not exceed allowed size for integers
        *lp->gslFlag_losint -= 1;   }
    return intg;
    
}
