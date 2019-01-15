#ifndef SIM_STRUCTS_H_INCLUDED
#define SIM_STRUCTS_H_INCLUDED

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include "sim-header.h"

struct bbhutil_params {
  str filename = "file";
  VD *full_field = NULL;
  VD wr_field;
  char *file = NULL;
  int *shape = NULL;
  int rank = 1;
  dbl *coords = NULL;
  dbl *data = NULL;
  bool write_this = false;
} ;
typedef struct bbhutil_params BBHP;

struct sim_fields {
  VD Al;
  VD oldAl;
  VD cnAl;
  VD resAl;
  VD olderAl;
  VD Be;
  VD oldBe;
  VD cnBe;
  VD resBe;
  VD olderBe; 
  VD Ps;
  VD oldPs;
  VD cnPs;
  VD resPs;
  VD olderPs;
  VD Xi;
  VD oldXi;
  VD cnXi;
  VD resXi;
  VD olderXi;
  VD Pi;
  VD oldPi;
  VD cnPi;
  VD resPi;
  VD olderPi;
  //VD res_hyp(n_hyp*npts, 0);
  VD res_ell;
  VD jac;
} ;
typedef struct sim_fields FLDS;

struct sim_writers {
  BBHP p_Al;
  BBHP p_resAl;
  BBHP p_iresAl;
  BBHP p_Be;
  BBHP p_resBe;
  BBHP p_iresBe;
  BBHP p_Ps;
  BBHP p_resPs;
  BBHP p_iresPs;
  BBHP p_Xi;
  BBHP p_resXi;
  BBHP p_iresXi;
  BBHP p_Pi;
  BBHP p_resPi;
  BBHP p_iresPi;
  BBHP p_maspect;
  BBHP p_outnull;
  BBHP p_ricci;
} ;
typedef struct sim_writers WRS;

struct sim_params {
  str outfile = "ekg";
  int lastpt = 500;
  int save_pt = 1;
  int nsteps = 1000;
  int save_step = 4;
  int resn_factor = 1;
  int check_step = 10;
  int maxit = 100;
  int norm_type = 0;
  dbl lam = 0.25; // dt/dr
  dbl rmin = 0;
  dbl rmax = 50;
  dbl r2m = 0;
  dbl dspn = 0.5; // dissipation coefficient
  dbl tol = 0.000000001; // iterative method tolerance
  dbl ell_tol = 0.01*tol;
  dbl ell_up_weight = 0.5;
  dbl ic_Dsq = 25.0; // gaussian width
  dbl ic_r0 = 10.0; // gaussian center
  dbl ic_Amp = 0.004; // gaussian amplitude
  bool psi_hyp = true; // update psi with hyperbolic evolution eqn after IC?
  bool zero_pi = false; // zero initial time derivative?
  bool somm_cond = true; // sommerfeld condition at outer bound?
  bool dspn_bound = false; // dissipate boundary points?
  bool dspn_psi = false; // dissipate psi (only activated if psi_hyp=true)?
  bool dr3_up = false; // update pi with d/dr^3 scheme?
  bool static_metric = false; // ignore scalar field's effect on metric?
  bool clean_hyp = false; // use clean hyperbolic update functions (slower)?
  bool clean_ell = false; // use clean hyperbolic update functions (slower)?
  bool write_res = false; // write residuals?
  bool write_ricci = false; // write ricci?
  bool write_itn = false; // write itn counts?
  bool write_mtot = false; // write total mass?
  bool write_maspect = false; // write mass aspect?
  bool write_outnull = false; // write outgoing null expansion?
  bool write_xp = false; // write xi and pi?
  bool write_abp = false; // write metric fields (alpha, beta, psi)?
  bool write_ires_xp = false; // write ires for xi and pi?
  bool write_ires_abp = false; // write ires for metric variables?
  bool horizon_search = true; // search for apparent horizon after each step?

  int n_ell = 3;
  int n_hyp = 2;

  int lastwr = lastpt/save_pt;
  int wr_shape = lastwr + 1;
  dbl wr_dr = (rmax - rmin) / ((dbl) lastwr);
  dbl coord_lims[2] = {rmin, rmax};
  vector< pair<int,int> > inds;
  int check_diagnostics = save_step * check_step;
  
  MAPID r {{0, rmin}};
  dbl t = 0;

  // DERIVED PARAMS
  int npts = lastpt + 1;
  dbl norm_factor = 1 / ((dbl) npts);
  dbl dr = (rmax - rmin) / ((dbl) lastpt);
  dbl dt = lam * dr;

  // LAPACK PARAMS
  int lp_n = n_ell * npts;
  int lp_kl = 2;
  int lp_ku = 2;
  int lp_nrhs = 1;
  int lp_ldab = 2*lp_kl + lp_ku + 1;
  int lp_ldb = lp_n;
  int *lp_ipiv = NULL;
  vector<int> ipiv; 
  
  // FREQUENTLY USED
  dbl one_third = 1.0 / 3.0;
  dbl two_thirds = 2 * one_third;
  dbl four_thirds = 2 * two_thirds;
  dbl twelfth = 0.25 * one_third;
  dbl five_twelfths = 5 * twelfth;
  dbl eight_pi = 8 * M_PI;
  dbl twelve_pi = 1.5 * eight_pi;
  dbl lam2val = 0.5 * lam;
  dbl lam6val = lam2val * one_third;
  dbl drsq = dr * dr;
  dbl indr = 1 / dr;
  dbl in2dr = 0.5 / dr;
  dbl indrsq = indr * indr;
  dbl neg2indrsq = -2 * indrsq;
  dbl indt = 1 / dt;
  dbl inrmax = 1 / rmax;
  // SPECIFIC TERMS
  dbl csomm = 0.75*lam + 0.5*dt*inrmax;
  dbl jacRR = 3*in2dr + inrmax;
  dbl jacRRm1 = -4 * in2dr;
  dbl jacRRm2 = in2dr;
  dbl jacN00 = -3 * in2dr;
  dbl jacN01 = 4 * in2dr;
  dbl jacN02 = -1 * in2dr;
  dbl dt_twelve = twelfth * dt;
  dbl cpsi_rhs = 1 / jacRR;
  dbl csomm_rhs = 1 / (1 + csomm);
  dbl csomm_old = 1 - csomm;
  
} ;
typedef struct sim_params PAR;

  
  
  





#endif
