/* ----------------------------------------------------------------------------
  CJAM
  
    For a given input model, the program reads in data-file locations and then
    calculates the moments of the model.  The moments are saved to specified
    file.
    
    INPUTS
      beta  : velocity anisotropy
      kappa : rotation parameter
      ml    : mass-to-light ratio [Msun/Lsun]
      incl  : inclination angle [radians]
      dist  : distance [kpc]
      mbh   : black-hole mass [Msun]
      rbh   : black-hole scale length [arcsec]
      nlg   : number of components in luminous MGE
      nmg   : number of components in mass MGE
      flmge : path to luminous MGE file
      fmmge : path to mass MGE file
      fxy   : path to star positions file
      fmom  : path to output moments file
  
  Laura L Watkins [lauralwatkins@gmail.com]
  
  This code is released under a BSD 2-clause license.
  If you use this code for your research, please cite:
  Watkins et al. 2013, MNRAS, 436, 2598
  "Discrete dynamical models of omega Centauri"
  http://adsabs.harvard.edu/abs/2013MNRAS.436.2598W
---------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "jam/jam.h"
#include "mge/mge.h"
#include "tools/readcol.h"
#include "tools/tools.h"


void cjam( double *beta, double *kappa, double *ml, double incl, double dist,
    double mbh, double rbh, int nlg, int nmg, char *flmge, char *fmmge,
    char *fxy, char *fmom ) {
    
    FILE *fp;
    struct multigaussexp lum, pot;
    char flags[60];
    int nrad, nang, i, j, k, filelen, nxy, check;
    double rad2arcsec, arcsec2pc, *xp, *yp;
    struct jam_vel vm;
    double *rxxm, *ryym, *rzzm, *rxym, *rxzm, *ryzm;
    
    
    
    // set flags for file reading
    filelen = 600000;
    sprintf( flags, " quiet bufsize=1024 filelen=%i skipsym=# ", filelen );
    
    // read in luminous MGE
    printf( "Reading luminous MGE from %s\n", flmge );
    mge_read( flmge, nlg, &lum );
    
    // read in potential MGE
    printf( "Reading potental MGE from %s\n", fmmge );
    mge_read( fmmge, nmg, &pot );
    
    // read in positions
    printf( "Reading positions from %s\n", fxy );
    xp = (double *) malloc( filelen * sizeof( double ) );
    yp = (double *) malloc( filelen * sizeof( double ) );
    nxy = readcol( fxy, flags, "%lf %lf", xp, yp );
    
    
    
    // convert arcsec to parsec
    rad2arcsec = 60. * 60. * RA2DEG;                    // radian -> arcsec
    arcsec2pc = dist * 1e3 / rad2arcsec;                // arcsec -> pc
    for ( i = 0; i < nxy; i++ ) {
        xp[i] *= arcsec2pc;
        yp[i] *= arcsec2pc;
    }
    for ( i = 0; i < lum.ntotal; i++ ) lum.sigma[i] *= arcsec2pc;
    for ( i = 0; i < pot.ntotal; i++ ) pot.sigma[i] *= arcsec2pc;
    
    // adjust potential mge by M/L
    for ( i = 0; i < pot.ntotal; i++ ) pot.area[i] *= ml[i];
    
    // add BH to potential gaussian
    rbh *= arcsec2pc;
    pot = mge_addbh( &pot, mbh, rbh );
    
    
    
    // properties of interpolation grid
    nrad = 25;
    nang = 10;
    
    // check for at least 1 rotating, non-spherical, non-isotropic component
    check = 0;
    for ( k = 0; k < lum.ntotal; k++ ) for ( j = 0; j < pot.ntotal; j++ ) 
        if ( ( kappa[k] != 0. ) & ( ( beta[k] != 0. ) | ( lum.q[k] != 1. ) 
            | ( pot.q[j] != 1. ) ) ) check++;
    
    if ( check > 0 ) {
        
        // calculate first moments
        printf( "Calculating first moments.\n" );
        vm = jam_axi_vel_mmt( xp, yp, nxy, incl, &lum, &pot, beta, kappa, nrad,
            nang );
    
    } else {
        printf( "Not calculating first moments -- model does not contain a "
            "rotating, non-spherical, non-isotropic component.\n" );
    }
    
    // calculate second moments
    printf( "Calculating second moments.\n" );
    rxxm = jam_axi_rms_mmt( xp, yp, nxy, incl, &lum, &pot, beta, \
        nrad, nang, 1 );
    ryym = jam_axi_rms_mmt( xp, yp, nxy, incl, &lum, &pot, beta, \
        nrad, nang, 2 );
    rzzm = jam_axi_rms_mmt( xp, yp, nxy, incl, &lum, &pot, beta, \
        nrad, nang, 3 );
    rxym = jam_axi_rms_mmt( xp, yp, nxy, incl, &lum, &pot, beta, \
        nrad, nang, 4 );
    rxzm = jam_axi_rms_mmt( xp, yp, nxy, incl, &lum, &pot, beta, \
        nrad, nang, 5 );
    ryzm = jam_axi_rms_mmt( xp, yp, nxy, incl, &lum, &pot, beta, \
        nrad, nang, 6 );
    
    
    
    // output model moments to file
    printf( "Writing moments to file %s\n\n", fmom );
    fp = fopen( fmom, "w" );
    for ( i = 0; i < nxy; i++ ) {
        if ( check == 0 ) fprintf( fp, "%lf  %lf  %lf  ", 0., 0., 0. );
        else fprintf( fp, "%lf  %lf  %lf  ", vm.vx[i], vm.vy[i], vm.vz[i] );
        fprintf( fp, "%lf  %lf  %lf  ", rxxm[i], ryym[i], rzzm[i] );
        fprintf( fp, "%lf  %lf  %lf\n", rxym[i], rxzm[i], ryzm[i] );
    }
    fclose( fp );
    
    
    
    // free memory
    free( xp );
    free( yp );
    free( lum.area );
    free( lum.sigma );
    free( lum.q );
    free( pot.area );
    free( pot.sigma );
    free( pot.q );
    
}
