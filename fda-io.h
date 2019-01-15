#ifndef FDA_IO_H_INCLUDED
#define FDA_IO_H_INCLUDED

#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <vector> // for everything
#include <cmath> // for ICs
#include <string> // for parameter input
#include <sstream> // for printing doubles to full precision
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf
#include "sim-structs.h"
#include "ekg-header.h"
#include "fda-fns.h"
#include "ekg-fns.h"
#include "ekg-proc.h"

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

inline void write_bbhp(BBHP *bp, PAR *p)
{
  prepare_write(*(bp->full_field), bp->wr_field, p);
  gft_out_bbox(bp->file, p->t, bp->shape, bp->rank, bp->coords, bp->data);
}

inline void write_bbhp_vec(vector<BBHP *>& bp_vec, PAR *p) {
  for (BBHP *bp : bp_vec) { write_bbhp(bp, p, p->t); }
}

VD make_vector(int len, dbl val) {
  VD new_vec(len, val);
  return new_vec;
}

void param_collect(char **source, int num, map<str, str>& dest) {
  for (int arg = 1; arg < num; ++arg) {
    if (source[arg][0] == '-') {
      dest[source[arg]] = source[arg+1];
    }
  }
}

void param_set(map<str, str>& p_all, map<str, str *>& p_str,
	       map<str, int *>& p_int, map<str, dbl *>& p_dbl,
	       map<str, bool *>& p_bool) {
  for (pair<str, str> p : p_all) {
    if (p_str.count(p.first)) { *p_str[p.first] = p.second; }
    else if (p_int.count(p.first)) { *p_int[p.first] = atoi(&p.second[0]); }
    else if (p_dbl.count(p.first)) { *p_dbl[p.first] = atof(&p.second[0]); }
    else if (p_bool.count(p.first)) { *p_bool[p.first] = (bool) atoi(&p.second[0]); }
  }
}

// read fields using bbhutil
void read_step(const vector<char *>& files, int times[], const vector<dbl *>& fields, int nfields) {
  for (int k = 0; k < nfields; ++k) {
    gft_read_brief(files[k], times[k], fields[k]);
  }
  return;
}

void write_diagnostics(WRS *wr, FLDS *f, PAR *p)
{
  // write ires
  // write ricci
  if (p->write_ricci) {
    get_ricci((wr->ricci).wr_field, f->Xi, f->Pi, f->Ps, p->inds);
    write_sdf(&(wr->ricci), p->t);
  }
  // write outnull
  if (p->write_outnull) {
    get_outnull((wr->outnull).wr_field, f->Al, f->Be, f->Ps,
		p->r, p->lastwr, p->save_pt);
    write_sdf(&(wr->outnull), p->t);
  }
  // write maspect
  if (p->write_maspect) {
    get_maspect((wr->maspect).wr_field, f->Al, f->Be, f->Ps,
		p->r, p->lastwr, p->save_pt);
    write_sdf(&(wr->maspect), p->t);
  }
  return;
}

/*
// write fields using bbhutil
void wr_step(int nfields, const vector<char *>& files,
	     dbl time, int *shape, int rank, dbl *coordinates,
	     const vector<dbl *>& fields) {
  for (int k = 0; k < nfields; ++k) {
    gft_out_bbox(files[k], time, shape, rank, coordinates, fields[k]);
  }
  return;
}

// GET COARSENED ARRAY FOR WRITING
void get_wr_f(const VD& f, VD& wr, int wr_shape, int save_pt)
{
  int k, s = 0;
  for (k = 0; k < wr_shape; ++k) {
    wr[k] = f[s];
    s += save_pt;
  }
  return;
}
*/

void record_horizon(PAR *p, const VD& f_ps, int ind, int itn, int t_itn)
{
  cout << "\nHORIZON FOUND\n" << endl;
  str param_str = "\noutfile name = " + p->outfile + "\nlaspt (save_pt) = " +
    to_string(p->lastpt) + " (" + to_string(p->save_pt) + ")\nnsteps (save_step) = " +
    to_string(p->nsteps) + " (" + to_string(p->save_step) + ")\nrmin = " + to_string(p->rmin)
    + ";\trmax = " + to_string(p->rmax) + "\nlambda = " + to_string(p->lam) +
    "\ndt = " + to_string(p->dt) + "\ndissipation = " + to_string(p->dspn) +
    "\nell_up_weight = " + to_string(p->ell_up_weight) + "\nhyp_tol = " + to_string(p->tol) +
    "\nell_tol = " + to_string(p->ell_tol) + "\nmaxit = " + to_string(maxit) + "\nic_r0 = " +
    to_string(p->ic_r0) + "\nic_Amp = " + to_string(p->ic_Amp) + "\nic_Dsq = " +
    to_string(p->ic_Dsq) + "\ndr = " + to_string(p->dr);
  param_str += "\noptions:\nhyperbolic psi evolution = " + bool_to_str(p->psi_hyp) +
    "\nzero pi_0 = " + bool_to_str(p->zero_pi) + "\nsommerfeld bc = " + bool_to_str(p->somm_cond)
    + "\ndissipation at bound = " + bool_to_str(p->dspn_bound) + "\nclean hyperbolic update functions = "
    + bool_to_str(p->clean_hyp) + "\nclean elliptic update functions = " + bool_to_str(p->clean_ell) + "\n";
  ofstream specs;
  str specs_name = "horizon-" + p->outfile + ".txt";
  specs.open(specs_name, ofstream::out);
  specs << "Horizon Found at:\nr[" << ind << "] = " << p->r[ind] << "  (r_areal = "
	<< sq(f_ps[ind])*p->r[ind] << ")\nt[" << t_itn << "] = "
	<< p->t << "  (itn " << itn << ")\nUsing:\ndr = " << p->dr << "\nic_Amp = " << p->ic_Amp
	<< "\nic_Dsq = " << p->ic_Dsq << "\nic_r0 = " << p->ic_r0 << "\n\n\nFULL PARAMETER DATA:\n";
  specs << param_str;
  specs.close();
  return;
}

#endif


