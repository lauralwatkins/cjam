/* ----------------------------------------------------------------------------
  CJAM_MAIN
  
    Parses input from command line and passes it to cjam function.
    
    INPUTS
      flmge   : path to luminous MGE
      nlg     : number of luminous MGE components
      fmmge   : path to mass MGE
      nmg     : number of mass MGE components
      fxy     : path to star positions file
      fmom    : path to ouput moments file
      incl    : inclination angle [radians]
      dist    : distance [kpc]
      mbh     : black-hole mass [Msun]
      rbh     : black-hole scale length [arcsec]
      (beta)  : velocity anisotropy (for each of the nlg components)
      (kappa) : rotation parameter (for each of the nlg components)
      (ml)    : mass-to-light ratio (for each of the nmg components)
      verbose : print progress, if set
  
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


void cjam( double*, double*, double*, double, double, double, double, int, int,
    char*, char*, char*, char*, int );


int main( int argc, char *argv[] ) {
    
    char *flmge, *fmmge, *fxy, *fmom;
    double incl, dist, mbh, rbh, *beta, *kappa, *ml, temp;
    int nlg, nmg, i, count, narg, verbose;
    
    
    
    if ( argc == 1 ) {
        printf( "\nRequired arguments:\n"
            "  flmge   : path to luminous MGE\n"
            "  nlg     : number of luminous MGE components\n"
            "  fmmge   : path to mass MGE\n"
            "  nmg     : number of mass MGE components\n"
            "  fxy     : path to star positions file\n"
            "  fmom    : path to ouput moments file\n"
            "  incl    : inclination angle [radians]\n"
            "  dist    : distance [kpc]\n"
            "  mbh     : black-hole mass [Msun]\n"
            "  rbh     : black-hole scale length [arcsec]\n"
            "  *beta   : velocity anisotropy (nlg values)\n"
            "  *kappa  : rotation parameter (nlg values)\n"
            "  *ml     : mass-to-light ratio (nmg values)\n"
            "  verbose : print progress, if set\n"
            "\nNote that beta and kappa are required for each\ncomponent of "
            "the luminous MGE, and ml is required for\neach component of the "
            "potential MGE.  In total there\nshould be (10 + 2*nlg + nmg) "
            "arguments provided.\n\n" );
        return 0;
    }
    
    
    
    // read verbose parameter
    sscanf( argv[argc-1], "%i", &verbose );
    
    if ( verbose ) printf( "\nReading parameters from command line...\n" );
    
    // check we have at least 4 parameters given so we can read the MGEs
    if ( argc - 1 < 4 ) {
        printf( "\nERROR: You have not provided enough arguments to read in "
        "the MGE files.\n\n" );
        return 0;
    }
    
    // read path to luminous MGE and size of luminous MGE
    flmge = argv[1];
    sscanf( argv[2], "%i", &nlg );
    if ( verbose ) printf( "  luminous MGE:    %s (%i components)\n", flmge, nlg );
    
    // read path to mass MGE and size of mass MGE
    fmmge = argv[3];
    sscanf( argv[4], "%i", &nmg );
    if ( verbose ) printf( "  potential MGE:   %s (%i components)\n", fmmge, nmg );
    
    
    
    // use nlg and nmg to check that correct number of betas/kappas/mls given
    narg = 11 + nlg * 2 + nmg;
    if ( argc - 1 < narg ) {
        printf( "\nERROR: Not enough arguments.  Expected %i, got %i.\n\n",
            narg, argc - 1 );
        return 0;
    }
    if ( argc - 1 > narg ) {
        printf( "\nERROR: Too many arguments.  Expected %i, got %i.\n\n",
            narg, argc - 1 );
        return 0;
    }
    
    
    
    // read path to positions file
    fxy = argv[5];
    if ( verbose ) printf( "  input positions: %s\n", fxy );
    
    // read path to moments file
    fmom = argv[6];
    if ( verbose ) printf( "  output moments:  %s\n", fmom );
    
    
    
    // read inclination, distance, black hole mass and black hole scale length
    sscanf( argv[7], "%lf", &incl );
    if ( verbose ) printf( "  inclination:     %lf radians\n", incl );
    sscanf( argv[8], "%lf", &dist );
    if ( verbose ) printf( "  distance:        %lf kpc\n", dist );
    sscanf( argv[9], "%lf", &mbh );
    if ( verbose ) printf( "  BH mass:         %lf Msun\n", mbh );
    sscanf( argv[10], "%lf", &rbh );
    if ( verbose ) printf( "  BH scale length: %lf arcsec\n", rbh );
    
    
    
    // counter to track arguments (because number of betas/kappas/mls can vary)
    count = 11;
    
    // read beta values into an array of size nlg
    beta = (double *) malloc( nlg * sizeof( double ) );
    if ( verbose ) printf( "  anisotropy:    " );
    for ( i = 0; i < nlg; i++ ) {
        sscanf( argv[count], "%lf", &temp );
        beta[i] = temp;
        if ( verbose ) printf( "  %lf", beta[i] );
        count++;
    }
    
    // read kappa values into an array of size nlg
    kappa = (double *) malloc( nlg * sizeof( double ) );
    if ( verbose ) printf( "\n  rotation:      " );
    for ( i = 0; i < nlg; i++ ) {
        sscanf( argv[count], "%lf", &temp );
        kappa[i] = temp;
        if ( verbose ) printf( "  %lf", kappa[i] );
        count++;
    }
    
    // read mass-to-light values into an array of size nmg
    ml = (double *) malloc( nmg * sizeof( double ) );
    if ( verbose ) printf( "\n  mass-to-light: " );
    for ( i = 0; i < nmg; i++ ) {
        sscanf( argv[count], "%lf", &temp );
        ml[i] = temp;
        if ( verbose ) printf( "  %lf", ml[i] );
        count++;
    }
    if ( verbose ) printf( "\n\n" );
    
    
    
    // call CJAM function
    cjam( beta, kappa, ml, incl, dist, mbh, rbh, nlg, nmg, flmge, fmmge, fxy,
        fmom, verbose );
    
    
    return 0;
}
