/* ----------------------------------------------------------------------------
  JAM_AXI_RMS_MGEINT
    
    Integrand for the MGE integral required for second moment calculation.
    
    INPUTS
      u      : integration variable
      params : function parameters passed as a structure
    
    NOTES
      * Based on janis2_jeans_mge_integrand IDL code by Michele Cappellari.
    
  Mark den Brok
  Laura L Watkins [lauralwatkins@gmail.com]
  
  This code is released under a BSD 2-clause license.
  If you use this code for your research, please cite:
  Watkins et al. 2013, MNRAS
  "Discrete dynamical models of omega Centauri"
  http://adsabs.harvard.edu/abs/2013MNRAS.tmp.2480W
---------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../mge/mge.h"
#include "jam.h"


double jam_axi_rms_mgeint( double u, void *params ) {
    
    struct params_rmsint *p;
    double e2u2p, a, b, c, d, e, f, sum;
    int j, k;
    
    double u2 = u * u;
    
    p = params;
    
    sum = 0.;
    for ( j = 0; j < p->pot->ntotal; j++ ) { //mass gaussians
        
        e2u2p = u2 * p->e2p[j];
        
        for ( k = 0; k < p->lum->ntotal; k++ ) { // luminous gaussians
            
            a = 0.5 * ( u2 / p->s2p[j] + 1. / p->s2l[k] );
            b = 0.5 * ( e2u2p * u2 / ( p->s2p[j] * ( 1. - e2u2p ) ) \
                + ( 1. - p->q2l[k] ) / p->s2q2l[k] );
            c = p->e2p[j] - p->s2q2l[k] / p->s2p[j];
            d = 1. - p->kani[k] * p->q2l[k] \
                - ( ( 1. - p->kani[k] ) * c + p->e2p[j] * p->kani[k] ) * u2;
            e = a + b * p->ci2;
            
            switch ( p->vv ) {
                case 1: // v2xx
                    f = p->kani[k] * p->s2q2l[k] + 0.5 * d * p->si2 / e \
                        + d * pow( a + b, 2 ) * p->ci2 * p->y2 / e / e;
                    break;
                case 2: // v2yy
                    f = p->s2q2l[k] * ( p->si2 + p->kani[k] * p->ci2 ) \
                        + p->x2 * p->ci2 * d;
                    break;
                case 3: // v2zz
                    f = p->s2q2l[k] * ( p->ci2 + p->kani[k] * p->si2 ) \
                        + p->x2 * p->si2 * d;
                    break;
                case 4: // v2xy
                    f = d * fabs( p->xy ) * p->ci2 * ( a + b ) / e;
                    break;
                case 5: // v2xz
                    f = d * fabs( p->xy ) * p->cisi * ( a + b ) / e;
                    break;
                case 6: // v2yz
                    f = p->cisi * ( p->s2q2l[k] * ( 1 - p->kani[k] ) \
                        - d * p->x2 );
                    break;
                default:
                    printf( "No integral selected.  Options: 1=v2xx, " );
                    printf( "2=v2yy, 3=v2zz, 4=v2xy, 5=v2xz, 6=v2yz.\n" );
                    f = 1.;
                    break;
            }
            
            sum += p->lum->area[k] * p->pot->q[j] * p->pot->area[j] * u2 \
                / ( 1. - c * u2 ) / sqrt( ( 1. - e2u2p ) * e ) \
                * f * exp( -a * ( p->x2 + p->y2 * ( a + b ) / e ) );
            
        }
    }
    
    return 4. * pow( M_PI, 1.5 ) * G * sum;
    
}
