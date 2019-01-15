#ifndef SIM_INIT_H_INCLUDED
#define SIM_INIT_H_INCLUDED

#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <cmath> // for ICs
#include <vector> // for everything
#include <string> // for parameter input
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf
#include "lapacke.h"
#include "fda-io.h"
#include "fda-fns.h"
#include "ekg-fns.h"
#include "jacobian.h"
//#include "ekg-clean.h"
#include "ekg-proc.h"
#include "solvers.h"
#include "sim-header.h"
#include "sim-structs.h"

void bbhp_init(BBHP *bp, PAR *p, str fieldname, VD *p_field, const VD& zeros)
{
  bp->full_field = p_field;
  bp->wr_field = zeros;
  str resn_str = "-" + to_string(p->resn_factor) + "-";
  bp->filename = fieldname + resn_str + p->outfile + ".sdf";
  bp->file = &(bp->filename[0]);
  bp->shape = &(p->wr_shape);
  bp->rank = 1;
  bp->coords = &(p->coord_lims[0]);  
  bp->data = &(bp->wr_field[0]);
  bp->write_this = true;
  return;
}

vector<BBHP *> writers_init(WRS *wr, FLDS *f, PAR *p)
{
  VD zeros(p->wr_shape, 0);
  vector<BBHP *> out_vec;
  if (p->write_abp) {
    bbhp_init(&(wr->p_Al), p, "Al", &(f->Al), zeros);
    bbhp_init(&(wr->p_Be), p, "Be", &(f->Be), zeros);
    bbhp_init(&(wr->p_Ps), p, "Ps", &(f->Ps), zeros);
    out_vec.push_back(&(wr->p_Al));
    out_vec.push_back(&(wr->p_Be));
    out_vec.push_back(&(wr->p_Ps));
    if (p->write_res) {
      bbhp_init(&(wr->p_resAl), p, "ResAl", &(f->resAl), zeros);
      bbhp_init(&(wr->p_resBe), p, "ResBe", &(f->resBe), zeros);
      bbhp_init(&(wr->p_resPs), p, "ResPs", &(f->resPs), zeros);
      out_vec.push_back(&(wr->p_resAl));
      out_vec.push_back(&(wr->p_resBe));
      out_vec.push_back(&(wr->p_resPs));
    }
  }

  if (p->write_xp) {
    bbhp_init(&(wr->p_Xi), p, "Xi", &(f->Xi), zeros);
    bbhp_init(&(wr->p_Pi), p, "Pi", &(f->Pi), zeros);
    out_vec.push_back(&(wr->p_Xi));
    out_vec.push_back(&(wr->p_Pi));
    if (p->write_res) {
      bbhp_init(&(wr->p_resXi), p, "ResXi", &(f->resXi), zeros);
      bbhp_init(&(wr->p_resPi), p, "ResPi", &(f->resPi), zeros);
      out_vec.push_back(&(wr->p_resXi));
      out_vec.push_back(&(wr->p_resPi));
    }
  }

  if (p->write_ires_abp) {
    bbhp_init(&(wr->p_iresAl), p, "iresAl", NULL, zeros);
    bbhp_init(&(wr->p_iresBe), p, "iresBe", NULL, zeros);
    bbhp_init(&(wr->p_iresPs), p, "iresPs", NULL, zeros);
  }
  if (p->write_ires_xp) {
    bbhp_init(&(wr->p_iresXi), p, "iresXi", NULL, zeros);
    bbhp_init(&(wr->p_iresPi), p, "iresPi", NULL, zeros);
  }
  if (p->write_maspect) { bbhp_init(&(wr->p_maspect), p, "maspect", NULL, zeros); }
  if (p->write_outnull) { bbhp_init(&(wr->p_outnull), p, "outnull", NULL, zeros); }
  if (p->write_ricci) { bbhp_init(&(wr->p_ricci), p, "ricci", NULL, zeros); }
  return out_vec;
}

int fields_init(FLDS *f, PAR *p)
{
  VD zeros(p->npts, 0);
  f->Al = zeros;
  f->Be = zeros;
  f->Ps = zeros;
  f->Xi = zeros;
  f->Pi = zeros;
  // PUT INITIAL CONDITIONS
  for (int k = 0; k < p->npts; ++k) {
    f->Al[k] = 1;
    f->Be[k] = 0;
    f->Ps[k] = 1;
    f->Xi[k] = ic_xi(p->r[k], p->ic_Amp, p->ic_Dsq, p->ic_r0);
    if (!(p->zero_pi)) { f->Pi[k] = ic_pi(p->r[k], p->ic_Amp, p->ic_Dsq, p->ic_r0); }
  }
  if (!(p->static_metric)) {
    int itn = solve_t0_slow(f->Xi, f->Pi, f->Al, f->Be, f->Ps, p, p->lastpt);
    if (itn < 0) {
      if (itn > -(p->npts)) {
	record_horizon(p, f->Ps, -itn, 0, 0);
      }
      else if (itn == -(p->npts)) {
	record_horizon(p, f->Ps, 0, 0, 0);
      }
      else { cout << "\nUNDETERMINED ERROR IN T = 0 ELLIPTIC CONSTRAINTS" << endl; }
      return itn;
    }
  }
  f->oldAl = f->Al;
  f->cnAl = f->Al;
  f->resAl = zeros;
  f->oldBe = f->Be;
  f->cnBe = f->Be;
  f->resBe = zeros;
  f->oldPs = f->Ps;
  f->cnPs = f->Ps;
  f->resPs = zeros;
  
  f->oldXi = f->Xi;
  f->cnXi = f->Xi;
  f->resXi = zeros;
  f->oldPi = f->Pi;
  f->cnPi = f->Pi;
  f->resPi = zeros;
  
  if (p->write_ires_abp) {
    f->olderAl = zeros;
    f->olderBe = zeros;
    f->olderPs = zeros;
  }
  if (p->write_ires_xp) {
    f->olderXi = zeros;
    f->olderPi = zeros;
  }

  VD res_zeros(p->lp_ldb, 0);
  VD jac_zeros(p->lp_ldab * p->lp_n, 0);
  // GET THESE ************************************
  f->res_ell = res_zeros;
  f->jac = jac_zeros;
  return 0;
}

int params_init(PAR *p, int argc, char **argv)
{
  map<str, str *> p_str {{"-outfile",&(p->outfile)}};
  map<str, int *> p_int {{"-lastpt",&(p->lastpt)}, {"-save_pt",&(p->save_pt)},
      {"-nsteps",&(p->nsteps)}, {"-save_step",&(p->save_step)},
      {"-norm_type",&(p->norm_type)}, {"-maxit",&(p->maxit)},
      {"-check_step",&(p->check_step)}, {"-resn_factor",&(p->resn_factor)}};
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
  if (p->lastpt % p->save_pt != 0) {
    cout << "ERROR: save_pt = " << p->save_pt << " entered for grid size " << p->lastpt << endl;
    p->save_pt -= p->lastpt % p->save_pt;
    cout << "--> corrected: save_pt = " << p->save_pt << endl;
  }
  // check that p->norm_type is valid
  if (p->norm_type < 0 || p->norm_type > 2) {
    cout << "ERROR: norm_type=" << p->save_pt << " not valid -> using inf-norm" << endl;
    p->norm_type = 0;
  }
  if (p->horizon_search) { cout << "\nsearching for horizon..." << endl; }
  //************************************************************************
  //********************SETTING PARAMS FROM OLD PROGRAM*********************
  if (p->psi_hyp) {
    p->n_ell = 2;
    p->n_hyp = 3;
  }
  // bbhutil parameters for writing data to sdf
  p->lastwr = p->lastpt / p->save_pt;
  p->wr_shape = p->lastwr + 1;
  p->coord_lims[0] = p->rmin; p->coord_lims[1] = p->rmax;
  p->wr_dr = (p->rmax - p->rmin) / ((dbl) p->lastwr);
  p->lastpt = p->lastpt * p->resn_factor;
  p->save_pt = p->save_pt * p->resn_factor;
  p->nsteps = p->nsteps * p->resn_factor;
  p->save_step = p->save_step * p->resn_factor;
  p->outfile = to_string(p->resn_factor) + "-" + p->outfile;
  // derived parameters
  p->npts = p->lastpt + 1;
  p->norm_factor = 1 / ((dbl) p->npts);
  if (p->norm_type == 1) { p->norm_factor = sqrt(p->norm_factor); }
  p->dr = (p->rmax - p->rmin) / ((dbl) p->lastpt);
  p->dt = p->lam * p->dr;
  p->check_diagnostics = p->save_step * p->check_step;
  
  for (int k = 0; k < p->npts; ++k) {
    p->r[k] = p->rmin + k * p->dr;
    if (p->r[k] != 0) { p->r[-k] = 1 / p->r[k]; }
  }
  p->t = 0;
  for (int k = 0; k < p->wr_shape; ++k) {
    (p->inds).push_back({k, (p->save_pt)*k});
  }
  if ((p->inds[p->lastwr]).second != p->lastpt) { cout << "\n***INDEX INIT ERROR***\n" << endl; }
  
  // lapack object declaration
  p->lp_n = p->n_ell * p->npts;
  p->lp_kl = 2;
  p->lp_ku = 2;
  p->lp_nrhs = 1;
  p->lp_ldab = 2 * p->lp_kl + p->lp_ku + 1;
  p->lp_ldb = p->lp_n;
  vector<lapack_int> ipiv_zeros(p->lp_n, 0);
  p->ipiv = ipiv_zeros;
  p->lp_ipiv = &(p->ipiv[0]);

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
    param_data += param.first + "\t\t=   " + bool_to_str(*(param.second)) + "\n";
  }
  ofstream specs;
  str specs_name = p->outfile + ".txt";
  specs.open(specs_name, ofstream::out);
  specs << param_data;
  specs.close();
  cout << param_data << endl;
  return 0;
}


#endif

