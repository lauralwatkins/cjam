/* ----------------------------------------------------------------------------
  JAM_AXI_VEL
    
    Wrapper for first moment calculator designed to interface with Python.
    
    INPUTS
      xp : projected x' [pc]
      yp : projected y' [pc]
      nxy : number of x' and y' values given
      incl : inclination [radians]
      lum_area : projected luminous MGE area
      lum_sigma : projected luminous MGE width
      lum_q : projected luminous MGE flattening
      pot_area : projected potential MGE area
      pot_sigma : projected potential MGE sigma
      pot_q : projected potential MGE flattening
      beta : velocity anisotropy (1-vz^2/vr^2)
      kappa : rotation parameter
      nrad : number of radial bins in interpolation grid
      nang : number of angular bins in interpolation grid
      vx : array to hold the vx first moments calculated
      vy : array to hold the vy first moments calculated
      vz : array to hold the vz first moments calculated
---------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "jam.h"
#include "../mge/mge.h"


void jam_axi_vel(double *xp, double *yp, int nxy, double incl, \
double *lum_area, double *lum_sigma, double *lum_q, int lum_total, \
double *pot_area, double *pot_sigma, double *pot_q, int pot_total, \
double *beta, double *kappa, int nrad, int nang, int* integrationFlag, \
double *vx, double *vy, double *vz) {
    
    struct multigaussexp lum, pot;
    struct jam_vel vm;
    int i, j, k, check;
    
    // put luminous MGE components into structure
    lum.area = lum_area;
    lum.sigma = lum_sigma;
    lum.q = lum_q;
    lum.ntotal = lum_total;
    
    // put potential MGE components into structure
    pot.area = pot_area;
    pot.sigma = pot_sigma;
    pot.q = pot_q;
    pot.ntotal = pot_total;
    
    // check for at least 1 rotating, non-spherical, non-isotropic component
    check = 0;
    for (k=0; k<lum.ntotal; k++) for (j=0; j<pot.ntotal; j++)
        if ( (kappa[k]!=0.) & ((beta[k]!=0.)|(lum.q[k]!=1.)|(pot.q[j]!=1.)) )
            check++;
    
    if (check>0) {
        // calculate moments and put into results arrays
        vm = jam_axi_vel_mmt(xp, yp, nxy, incl, &lum, &pot, beta, kappa,
            nrad, nang, integrationFlag);
        for (i=0; i<nxy; i++) {
            vx[i] = vm.vx[i];
            vy[i] = vm.vy[i];
            vz[i] = vm.vz[i];
        }
        
        // free memory
        free(vm.vx);
        free(vm.vy);
        free(vm.vz);
    } else {
        // return zeros
        for (i=0; i<nxy; i++) {
            vx[i] = 0.;
            vy[i] = 0.;
            vz[i] = 0.;
        }
    }
    
    return;
}
