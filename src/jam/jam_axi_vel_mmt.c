/* -----------------------------------------------------------------------------
  JAM_AXI_VEL_MMT
    
    Calculates first moments.
    
    INPUTS
      xp    : projected x' [pc]
      yp    : projected y' [pc]
      nxy   : number of x' and y' values given
      incl  : inclination [radians]
      lum   : projected luminous MGE
      pot   : projected potential MGE
      beta  : velocity anisotropy (1 - vz^2 / vr^2)
      kappa : rotation parameter
      nrad  : number of radial bins in interpolation grid
      nang  : number of angular bins in interpolation grid
    
    NOTES
      * Based on janis1_first_moment IDL code by Michele Cappellari.
      * This version does not implement PSF convolution.
    
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
#include "jam.h"
#include "../mge/mge.h"
#include "../tools/tools.h"
#include "../interp/interp.h"


struct jam_vel jam_axi_vel_mmt( double *xp, double *yp, int nxy, \
        double incl, struct multigaussexp *lum, struct multigaussexp *pot, \
        double *beta, double *kappa, int nrad, int nang ) {
    
    int i, j, k, v, npol;
    double qmed, *rell, *r, *e, step, rmax, *lograd, *rad, *ang, *angvec;
    double **wm1, *surf, *xpol, *ypol, **mupol, *temp;
    struct jam_vel mu;
    
    
    mu.vx = (double *) malloc( nxy * sizeof( double ) );
    mu.vy = (double *) malloc( nxy * sizeof( double ) );
    mu.vz = (double *) malloc( nxy * sizeof( double ) );
    
    // skip the interpolation when computing just a few points
    if ( nrad * nang > nxy ) {
        
        // weighted first moments
        wm1 = jam_axi_vel_wmmt( xp, yp, nxy, incl, lum, pot, beta, kappa );
        
        // surface brightness
        surf = mge_surf( lum, xp, yp, nxy );
        
        // first moments
        for ( i = 0; i < nxy; i++ ) {
            mu.vx[i] = wm1[i][0] / surf[i];
            mu.vy[i] = wm1[i][1] / surf[i];
            mu.vz[i] = wm1[i][2] / surf[i];
        }
        
        for ( i = 0; i < 3; i++ ) free( wm1[i] );
        free( wm1 );
        free( surf );
        
        return mu;
    }
    
    
    // ---------------------------------
    
    
    // elliptical radius of input (x,y)
    qmed = mge_qmed( lum, maximum( xp, nxy ) );
    rell = (double *) malloc( nxy * sizeof( double ) );
    for ( i = 0; i < nxy; i++ ) \
        rell[i] = sqrt( xp[i] * xp[i] + yp[i] * yp[i] / qmed / qmed );
    
    // set interpolation grid parameters
    step = minimum( rell, nxy );
    if ( step <= 0.001 ) step = 0.001;          // minimum radius of 0.001 pc
    npol = nrad * nang;
    
    // make linear grid in log of elliptical radius
    rmax = maximum( rell, nxy );
    lograd = range( log( step ) - 0.1, log( rmax ) + 0.1, nrad, False ); 
    rad = (double *) malloc( nrad * sizeof( double ) );
    for ( i = 0; i < nrad; i++ ) rad[i] = exp( lograd[i] );
    
    // make linear grid in eccentric anomaly
    ang = range( -M_PI, -M_PI / 2., nang, False );
    angvec = range( -M_PI, M_PI, 4 * nang - 3, False );
    
    // convert grid to cartesians
    xpol = (double *) malloc( npol * sizeof( double ) );
    ypol = (double *) malloc( npol * sizeof( double ) );
    for ( i = 0; i < nrad; i++ ) {
        for ( j = 0; j < nang; j++ ) {
            xpol[i*nang+j] = rad[i] * cos( ang[j] );
            ypol[i*nang+j] = rad[i] * sin( ang[j] ) * qmed;
        }
    }
    
    // set up interpolation grid arrays
    mupol = (double **) malloc( nrad * sizeof( double * ) );
    for ( i = 0; i < nrad; i++ ) \
        mupol[i] = (double *) malloc( ( 4 * nang - 3 ) * sizeof( double ) );
    
    
    // ---------------------------------
    
    
    // elliptical radius and eccentric anomaly of inputs
    r = (double *) malloc( nxy * sizeof( double ) );
    e = (double *) malloc( nxy * sizeof( double ) );
    for ( i = 0; i < nxy; i++ ) {
        r[i] = sqrt( pow( xp[i], 2. ) + pow( yp[i] / qmed, 2. ) );
        e[i] = atan2( yp[i] / qmed, xp[i] );
    }
    
    // weighted first moments on polar grid
    wm1 = jam_axi_vel_wmmt( xpol, ypol, npol, incl, lum, pot, beta, kappa );
    
    // surface brightness on polar grid
    surf = mge_surf( lum, xpol, ypol, npol );
    
    for ( v = 0; v < 3; v++ ) {
        
        // velocity first moment on polar grid
        for ( i = 0; i < nrad; i++ ) {
            for ( j = 0; j < nang; j++ ) {
                
                k = i * nang + j;
                mupol[i][j] = wm1[k][v] / surf[k];
                
                mupol[i][2*nang-2-j] = mupol[i][j];
                if ( v == 1 || v == 2 ) mupol[i][2*nang-2-j] *= -1.;
                mupol[i][2*nang-2+j] = -mupol[i][j];
                mupol[i][4*nang-4-j] = mupol[i][j];
                if ( v == 0 ) mupol[i][4*nang-4-j] *= -1;
                
            }
        }
        
        // interpolate to get first moment at input positions
        temp = interp2dpol( mupol, rad, angvec, r, e, nrad, 4*nang-3, nxy );
        
        if ( v == 0 ) mu.vx = temp;
        if ( v == 1 ) mu.vy = temp;
        if ( v == 2 ) mu.vz = temp;
        
    }
    
    
    // ---------------------------------
    
    
    free( rell );
    free( rad );
    free( ang );
    free( angvec );
    free( xpol );
    free( ypol );
    for ( i = 0; i < nrad; i++ ) free( mupol[i] );
    free( mupol );
    for ( i = 0; i < 3; i++ ) free( wm1[i] );
    free( wm1 );
    free( surf );
    free( r );
    free( e );
    
    return mu;
    
}
