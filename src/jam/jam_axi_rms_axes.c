/* ----------------------------------------------------------------------------
  JAM_AXI_RMS_AXES
    
    Wrapper for second moment calculator designed to interface with Python, calculating only requested moments.
    
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
      nrad : number of radial bins in interpolation grid
      nang : number of angular bins in interpolation grid
      rxx : array to hold the xx second moments calculated
      ryy : array to hold the yy second moments calculated
      rzz : array to hold the zz second moments calculated
      rxy : array to hold the xy second moments calculated
      rxz : array to hold the xz second moments calculated
      ryz : array to hold the yz second moments calculated
      xaxis : whether to calculate moments involving x
      yaxis : whether to calculate moments involving y
      zaxis : whether to calculate moments involving z
---------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "jam.h"
#include "../mge/mge.h"


void jam_axi_rms_axes(double *xp, double *yp, int nxy, double incl, \
double *lum_area, double *lum_sigma, double *lum_q, int lum_total, \
double *pot_area, double *pot_sigma, double *pot_q, int pot_total, \
double *beta, int nrad, int nang, int* integrationFlag, \
double *rxx, double *ryy, double *rzz, double *rxy, double *rxz, double *ryz,
int xaxis, int yaxis, int zaxis) {
    
    struct multigaussexp lum, pot;
    double* mu;
    int i, check;
    
    // if there are no moments requested, exit immediately
    if (!xaxis && !yaxis && !zaxis) {
        return;
    }
    
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
    
    // check for any non-zero beta or non-unity flattening
    check = 0;
    for (i=0; i<lum.ntotal; i++) if (beta[i]!=0.) check++;
    for (i=0; i<lum.ntotal; i++) if (lum_q[i]!=1.) check++;
    for (i=0; i<pot.ntotal; i++) if (pot_q[i]!=1.) check++;
    
    // for anisotropic models, calculate each requested moment
    if (check>0) {
        if (xaxis) {
            // calculate xx moments and put into results array
            mu = jam_axi_rms_mmt(xp, yp, nxy, incl, &lum, &pot, beta, nrad, nang, \
                1, integrationFlag);
            for (i=0; i<nxy; i++) rxx[i] = mu[i];
        }
        if (yaxis) {
            mu = jam_axi_rms_mmt(xp, yp, nxy, incl, &lum, &pot, beta, nrad, nang, \
                2, integrationFlag);
            for (i=0; i<nxy; i++) ryy[i] = mu[i];
        }
        if (zaxis) {
            mu = jam_axi_rms_mmt(xp, yp, nxy, incl, &lum, &pot, beta, nrad, nang, \
                3, integrationFlag);
            for (i=0; i<nxy; i++) rzz[i] = mu[i];
        }
        if (xaxis && yaxis) {
            mu = jam_axi_rms_mmt(xp, yp, nxy, incl, &lum, &pot, beta, nrad, nang, \
                4, integrationFlag);
            for (i=0; i<nxy; i++) rxy[i] = mu[i];
        }
        if (xaxis && zaxis) {
            mu = jam_axi_rms_mmt(xp, yp, nxy, incl, &lum, &pot, beta, nrad, nang, \
                5, integrationFlag);
            for (i=0; i<nxy; i++) rxz[i] = mu[i];
        }
        if (yaxis && zaxis) {
            mu = jam_axi_rms_mmt(xp, yp, nxy, incl, &lum, &pot, beta, nrad, nang, \
                6, integrationFlag);
            for (i=0; i<nxy; i++) ryz[i] = mu[i];
        }
    }
    // otherwise just calculate one and propagate
    else {
        mu = jam_axi_rms_mmt(xp, yp, nxy, incl, &lum, &pot, beta, nrad, nang, 1, \
            integrationFlag);
        if (xaxis) for (i=0; i<nxy; i++) rxx[i] = mu[i];
        if (yaxis) for (i=0; i<nxy; i++) ryy[i] = mu[i];
        if (zaxis) for (i=0; i<nxy; i++) rzz[i] = mu[i];
        if (xaxis && yaxis) for (i=0; i<nxy; i++) rxy[i] = 0.;
        if (xaxis && zaxis) for (i=0; i<nxy; i++) rxz[i] = 0.;
        if (yaxis && zaxis) for (i=0; i<nxy; i++) ryz[i] = 0.;
    }
    
    // free memory
    free(mu);
    
    return;
}
