#ifndef SIM_STRUCTS_H_INCLUDED
#define SIM_STRUCTS_H_INCLUDED

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include "sim-header.h"

struct bbhutil_params {
  MAPID *full_field;
  VD wr_field;
  str filename;
  char *file;
  int *shape;
  int rank;
  dbl *coords;
  dbl *data;
  bool write_this;
} ;
//typedef struct bbhutil_params BBHP;

struct sim_fields {
  MAPID Al;
  MAPID oldAl;
  MAPID resAl;
  MAPID olderAl;
  MAPID Be;
  MAPID oldBe;
  MAPID resBe;
  MAPID olderBe;  
  MAPID Ps;
  MAPID oldPs;
  MAPID resPs;
  MAPID olderPs;
  MAPID Xi;
  MAPID oldXi;
  MAPID resXi;
  MAPID olderXi;
  MAPID Pi;
  MAPID oldPi;
  MAPID resPi;
  MAPID olderPi;
  //VD res_hyp(n_hyp*npts, 0);
} ;
//typedef struct sim_fields FLDS;

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
//typedef struct sim_writers WRS;

struct sim_params {
  str outfile;
  int lastpt;
  int save_pt;
  int nsteps;
  int save_step;
  int resn_factor;
  int check_step;
  int maxit;
  int norm_type = 0;
  dbl lam = 0.25; // dt/dr
  dbl rmin = 0;
  dbl rmax;
  dbl r2m = 0;
  dbl dspn; // dissipation coefficient
  dbl tol; // iterative method tolerance
  dbl ell_tol;
  dbl ell_up_weight;
  dbl ic_Dsq = 25.0; // gaussian width
  dbl ic_r0 = 15.0; // gaussian center
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

  int n_ell;
  int n_hyp;

  int lastwr;
  int wr_shape;
  vector<int[2]> inds;
  dbl coor_lims[2];

  // DERIVED PARAMS
  int npts;
  dbl norm_factor;
  dbl dr;
  dbl dt;

  // LAPACK PARAMS
  int lp_n;
  int lp_kl;
  int lp_ku;
  int lp_nrhs = 1;
  int lp_ldab;
  int lp_ldb;
  int *ipiv;
  
  // FREQUENTLY USED
  dbl one_third = 1.0 / 3.0;
  dbl two_thirds = 2 * one_third;
  dbl four_thirds = 2 * two_thirds;
  dbl twelfth = 0.25 * one_third;
  dbl five_twelfths = 5 * twelfth;
  dbl eight_pi = 8 * M_PI;
  dbl twelve_pi = 1.5 * eight_pi;
  dbl lam2val;
  dbl lam6val;
  dbl indr;
  dbl in2dr;
  dbl drsq;
  dbl indrsq;
  dbl indt;
  dbl csomm;
  dbl jacRR;
  dbl jacRRm1;
  dbl jacRRm2;
  dbl jacN00;
  dbl jacN01;
  dbl jacN02;
  dbl dt_twelve;
  dbl cpsi_rhs;
  dbl csomm_rhd;
  dbl csomm_old;
  
} ;
//typedef struct sim_params PAR;

  
  
  





#endif
