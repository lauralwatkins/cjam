/* ----------------------------------------------------------------------------
  MGE_ADDBH
    
    Adds a black hole component to an MGE.
    
    INPUTS
      mge : MGE [area in M/pc^2, sigma in pc]
      mbh : black hole mass [Msun]
      rbh : black hole softening length [pc]
  
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
#include "mge.h"


struct multigaussexp mge_addbh( struct multigaussexp *mge, double mbh, \
    double rbh ) {
    
    struct multigaussexp mgebh;
    int i, isbh;
    
    // check if black hole needs to be added
    isbh = 0;
    if ( mbh > 0. && rbh > 0. ) isbh = 1;
    
    // memory allocation
    mgebh.ntotal = mge->ntotal + isbh;
    mgebh.sigma = ( double * ) malloc( mgebh.ntotal * sizeof( double ) );
    mgebh.area = ( double * ) malloc( mgebh.ntotal * sizeof( double ) );
    mgebh.q = ( double * ) malloc( mgebh.ntotal * sizeof( double ) );
    
    // add black hole MGE
    if ( isbh == 1 ) {
        mgebh.area[0] = mbh / 2. / M_PI / rbh / rbh;
        mgebh.sigma[0] = rbh;
        mgebh.q[0] = 1.;
    }
    
    // copy existing mass MGE
    for ( i = 0; i < mge->ntotal; i++ ) {
        mgebh.area[i+isbh] = mge->area[i];
        mgebh.sigma[i+isbh] = mge->sigma[i];
        mgebh.q[i+isbh] = mge->q[i];
    }
    
    return mgebh;
    
}
