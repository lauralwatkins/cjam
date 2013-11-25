/* -----------------------------------------------------------------------------
  JAM_AXI_VEL_MGEINT
    
    Calculates inner intergrand required for first moments.
    
    INPUTS
      u      : integration variable
      params : function parameters passed as a structure
      
    NOTES
      * Based on janis1_jeans_mge_integrand IDL code by Michele Cappellari.
    
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
#include "../mge/mge.h"
#include "jam.h"


double jam_axi_vel_mgeint( double u, void *params ) {
    
    struct params_mgeint *p;
    double c, d, p2, hj, e, sum;
    int j, k;
    
    double u2 = u * u;
    
    p = params;
    
    // double summation of eqn 38 over integration variable u
    sum = 0.;
    for ( j = 0; j < p->pot->ntotal; j++ ) { // mass gaussians
        
        p2 = 1. - p->e2p[j] * u2;
        hj = exp( -0.5 / p->s2p[j] * u2 * ( p->r2 + p->z2 / p2 ) ) \
            / sqrt( p2 );                                           // eqn 17
        e = p->pot->q[j] * p->pot->area[j] * hj * u2;
        
        for ( k = 0; k < p->lum->ntotal; k++ ) { // luminous gaussians
            
            c = p->e2p[j] - p->s2q2l[k] / p->s2p[j];                // eqn 22
            d = 1. - p->bani[k] * p->q2l[k] - ( ( 1. - p->bani[k] ) * c \
                + p->e2p[j] * p->bani[k] ) * u2;                    // eqn 23
            sum += p->knu[k] * e * d / ( 1. - c * u2 );             // eqn 38
            
        }
        
    }
    
    return sum;
    
}
