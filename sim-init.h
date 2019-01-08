#ifndef SIM_INIT_H_INCLUDED
#define SIM_INIT_H_INCLUDED

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <string>

#include "sim-header.h"
#include "sim-structs.h"
#include "fda-fns.h"
#include "fda-io.h"

inline str bool_to_str(bool bool_in)
{
  return ( (bool_in) ? "TRUE" : "FALSE" );
}


inline void write_sdf(const BBHP *bp, dbl t)
{
  gft_out_bbox(bp->file, t, bp->shape, bp->rank, bp->coords, bp->data);
}

inline void prepare_write(const MAPID& fld, VD& wr_fld, PAR *p)
{
  for (int ind[2] : p->inds) {
    wr_fld[ind[0]] = fld[ind[1]];
  }
}

inline void write_bbhp(BBHP *bp, PAR *p, dbl t)
{
  prepare_write(*(bp->full_field, bp->wr_field, p);
  gft_out_bbox(bp->file, t, bp->shape, bp->rank, bp->coords, bp->data);
}

void bbhp_init(BBHP *bp, PAR *p, str fieldname, MAPID *p_field, VD& zeros)
{
  bp->full_field = p_field;
  bp->wr_field = zeros;
  str resn_str = "-" + to_string(p->resn_factor) + "-";
  bp->filename = fieldname + resn_str + p->outfile + ".sdf";
  bp->file = &(bp->filename[0]);
  bp->shape = &(p->wr_shape);
  bp->rank = 1;
  bp->coords = &(p->coord_lims);  
  bp->data = &(bp->wr_field[0]);
  if (zeros.size() < 2) {
    bp->shape = &(bp->rank);
    bp->write_this = false;
  }
  else { bp->write_this = true; }
  return;
}

void writers_init(WRS *wr, FLDS *f, PAR *p)
{
  VD zeros(p->wr_shape, 0);
  VD unused(1, 0);
  bbhp_init(&(wr->p_Al), p, "Al", &(f->Al), ((p->write_abp) ? zeros : unused));
  bbhp_init(&(wr->p_Be), p, "Be", &(f->Be), ((p->write_abp) ? zeros : unused));
  bbhp_init(&(wr->p_Ps), p, "Ps", &(f->Ps), ((p->write_abp) ? zeros : unused));
  bbhp_init(&(wr->p_Xi), p, "Xi", &(f->Xi), ((p->write_xp) ? zeros : unused));
  bbhp_init(&(wr->p_Pi), p, "Pi", &(f->Pi), ((p->write_xp) ? zeros : unused));
  if (p->write_res) {
     bbhp_init(&(wr->p_resAl), p, "ResAl", &(f->resAl), ((p->write_abp) ? zeros : unused));
     bbhp_init(&(wr->p_resBe), p, "ResBe", &(f->resBe), ((p->write_abp) ? zeros : unused));
     bbhp_init(&(wr->p_resPs), p, "ResPs", &(f->resPs), ((p->write_abp) ? zeros : unused));
     bbhp_init(&(wr->p_resXi), p, "ResXi", &(f->resXi), ((p->write_xp) ? zeros : unused));
     bbhp_init(&(wr->p_resPi), p, "ResPi", &(f->resPi), ((p->write_xp) ? zeros : unused));
  }
  else {
    bbhp_init(&(wr->p_resAl), p, "ResAl", &(f->resAl), unused);
    bbhp_init(&(wr->p_resBe), p, "ResBe", &(f->resBe), unused);
    bbhp_init(&(wr->p_resPs), p, "ResPs", &(f->resPs), unused);
    bbhp_init(&(wr->p_resXi), p, "ResXi", &(f->resXi), unused);
    bbhp_init(&(wr->p_resPi), p, "ResPi", &(f->resPi), unused);
  }
  bbhp_init(&(wr->p_iresAl), p, "iresAl", &(f->iresAl), ((p->write_ires_abp) ? zeros : unused));
  bbhp_init(&(wr->p_iresBe), p, "iresBe", &(f->iresBe), ((p->write_ires_abp) ? zeros : unused));
  bbhp_init(&(wr->p_iresPs), p, "iresPs", &(f->iresPs), ((p->write_ires_abp) ? zeros : unused));
  bbhp_init(&(wr->p_iresXi), p, "iresXi", &(f->iresXi), ((p->write_ires_xp) ? zeros : unused));
  bbhp_init(&(wr->p_iresPi), p, "iresPi", &(f->iresPi), ((p->write_ires_xp) ? zeros : unused));

  bbhp_init(&(wr->p_maspect), p, "maspect", NULL, ((write_maspect) ? zeros : unused));
  bbhp_init(&(wr->p_outnull), p, "outnull", NULL, ((write_outnull) ? zeros : unused));
  bbhp_init(&(wr->p_ricci), p, "ricci", NULL, ((write_ricci) ? zeros : unused));
  return;
}

void fields_init(FLDS *f, PAR *p)
{
  // can incorporate options from p to save memory is, e.g., no ires
  VD zeros(p->npts, 0);
  VD unused(1, 0);
  f->Al = zeros;
  f->oldAl = zeros;
  f->resAl = zeros;
  f->olderAl = ((write_ires_abp) ? zeros : unused);

  f->Be = zeros;
  f->oldBe = zeros;
  f->resBe = zeros;
  f->olderBe = ((write_ires_abp) ? zeros : unused);
  
  f->Ps = zeros;
  f->oldPs = zeros;
  f->resPs = zeros;
  f->olderPs = ((write_ires_abp) ? zeros : unused);

  f->Xi = zeros;
  f->oldXi = zeros;
  f->resXi = zeros;
  f->olderXi = ((write_ires_xp) ? zeros : unused);

  f->Pi = zeros;
  f->oldPi = zeros;
  f->resPi = zeros;
  f->olderPi = ((write_ires_xp) ? zeros : unused);
  return;
}

str params_init(PAR *p, int argc, char **argv)
{
  map<str, str *> p_str {{"-outfile",&(p->outfile)}};
  map<str, int *> p_int {{"-lastpt",&(p->lastpt)}, {"-save_pt",&(p->save_pt)},
      {"-nsteps",&(p->nsteps)}, {"-save_step",&(p->save_step)},
      {"-norm_type",&(p->norm_type)}, {"-maxit",&(p->maxit)},
      {"-check_step",&(p->check_step)}, {"-resn_factor",&(p->nresn)}};
  map<str, dbl *> p_dbl {{"-lam",&(p->lam)}, {"-r2m",&(p->r2m)}, {"-rmin",&(p->rmin)}, {"-rmax",&(p->rmax)},
      {"-dspn",&(p->dspn)}, {"-tol",&(p->tol)}, {"-ell_tol",&(p->ell_tol)}, {"-ell_up_weight",&(p->ell_up_weight)},
      {"-ic_Dsq",&(p->ic_Dsq)}, {"-ic_r0",&(p->ic_r0)}, {"-ic_Amp",&(p->ic_Amp)}};
  map<str, bool *> p_bool { {"-psi_hyp",&(p->psi_hyp)}, {"-zero_pi",&(p->zero_pi)}, {"-static_metric",&(p->static_metric)},
      {"-somm_cond",&(p->somm_cond)}, {"-dspn_bound",&(p->dspn_bound)}, {"-dr3_up",&(p->dr3_up)}, {"-dspn_psi",&(p->dspn_psi)},
      {"-write_res",&(p->write_res)},{"-write_ricci",&(p->write_ricci)}, {"-write_itn",&(p->write_itn)},
      {"-write_mtot",&(p->write_mtot)},{"-write_maspect",&(p->write_maspect)}, {"-write_outnull",&(p->write_outnull)},
      {"-write_xp",&(p->write_xp)}, {"-write_abp",&(p->write_abp)},
      {"-write_ires_xp",&(p->write_ires_xp)}, {"-write_ires_abp",&(p->write_ires_abp)},
      {"-clean_hyp",&(p->clean_hyp)}, {"-clean_ell",&(p->clean_ell)}, {"-horizon_search",&(p->horizon_search)}};
  map<str, str> params;
  param_collect(argv, argc, params);
  param_set(params, p_str, p_int, p_dbl, p_bool);
  

  // *****************************************CHECKS**********************
  // check that grid size (lastpt = npts-1) is divisible by save_pt 
  if (lastpt % save_pt != 0) {
    cout << "ERROR: save_pt = " << save_pt << " entered for grid size " << lastpt << endl;
    save_pt -= lastpt % save_pt;
    cout << "--> corrected: save_pt = " << save_pt << endl;
  }
  // check that resnorm_type is valid
  if (resnorm_type < 0 || resnorm_type > 2) {
    cout << "ERROR: resnorm_type=" << save_pt << " not valid -> using inf-norm" << endl;
    resnorm_type = 0;
  }
  if (horizon_search) { cout << "\nsearching for horizon..." << endl; }
  //************************************************************************
  //********************SETTING PARAMS FROM OLD PROGRAM*********************
  int n_ell = 3;
  if (psi_hyp) { n_ell = 2; }
  int n_hyp = 5 - n_ell;
  // bbhutil parameters for writing data to sdf
  int lastwr = lastpt/save_pt;
  int wr_shape = lastwr + 1;
  int *bbh_shape = &wr_shape;
  int bbh_rank = 1;
  dbl coord_lims[2] = {rmin, rmax};
  dbl *coords = &coord_lims[0];
  dbl wr_dr = (rmax - rmin) / ((dbl) lastwr);
  lastpt = lastpt0 * factor;
  save_pt = save_pt0 * factor;
  nsteps = nsteps0 * factor;
  save_step = save_step0 * factor;
  outfile = to_string(factor) + "-" + outfile0;
  // derived parameters
  int npts = lastpt + 1;
  dbl norm_factor = 1 / ((dbl) 5*npts);
  if (resnorm_type == 1) { norm_factor = sqrt(norm_factor); }
  dbl dr = (rmax - rmin) / ((dbl) lastpt);
  dbl dt = lam * dr;
  MAPID r { {WRITEDR,wr_dr}, {HYPTOL,tol}, {ELLTOL,ell_tol},
					     {EUP_WEIGHT,ell_up_weight}, {DSPN_WEIGHT,dspn} };
  set_rmap(r, lastpt, dr, dt, lam, rmin, rmax);
  
  // lapack object declaration
  lapack_int N = n_ell*npts;
  lapack_int kl = 2;
  lapack_int ku = 2;
  lapack_int nrhs = 1;
  lapack_int ldab = 2*kl + ku + 1;
  lapack_int ldb = N;
  vector<lapack_int> ipiv(N);
  VD res_ell(ldb, 0);
  VD jac_zero(ldab*N, 0);
  
  MAPII z { {MAX_ITN,maxit}, {RESN_FACTOR,factor} };
  set_zmap(z, lastwr, lastpt, save_pt, nsteps, save_step, n_ell, kl, ku);
  //*****************************************************************************
  //*******************************************INDICES**************************
  for (int k = 0; k < p->wr_shape; ++k) {
    int ind[2] = {k, (p->save_pt)*k};
    (p->inds).push_back(ind);
  }
  if (inds[lastwr][1] != lastpt) { cout << "\n***INDEX INIT ERROR***\n" << endl; }


  // PARAMETER DATA OUTPUT
  str param_data = "\nPARAMETERS:\n\n";
  for (pair<str, str *> param : p_str) {
    param_data += param.first + "\t\t=   " + *(param.second) + "\n";
  }
  param_data += "\n";
  for (pair<str, int *> param : p_int) {
    param_data += param.first + "\t\t=   " + to_string(*(param.second)) + "\n";
  }
  param_data += "\n";
  for (pair<str, dbl *> param : p_dbl) {
    param_data += param.first + "\t\t=   " + to_string(*(param.second)) + "\n";
  }
  param_data += "\n";
  for (pair<str, bool *> param : p_bool) {
    param_data += param.first + "\t\t=   " + bool_to_string(*(param.second)) + "\n";
  }
  // THEN INIT THE INDICES VECTOR
  return param_data;

}


#endif

