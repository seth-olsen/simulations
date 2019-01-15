#ifndef SOLVERS_H_INCLUDED
#define SOLVERS_H_INCLUDED

#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <vector> // for everything
#include <cmath> // for ICs
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf
#include "lapacke.h"
#include "sim-header.h"
#include "sim-structs.h"
#include "fda-io.h"
#include "fda-fns.h"
#include "ekg-fns.h"
#include "jacobian.h"
#include "ekg-proc.h"
#include "ekg-clean.h"

using namespace std;

typedef int (*SOLVER)(VD& , VD& , VD& , VD& , VD& ,
		      VD& , VD& , VD& , VD& , VD& ,
		      VD& , VD& , VD& , VD& , VD& ,
		      VD& , VD& , const VD& , MAPID& ,
		      int , int , int ,
		      int , int , int , int , int , vector<int>& , int);

void fields_step(FLDS *f, PAR *p)
{
  f->oldAl = f->Al;
  f->cnAl = f->Al;  
  f->oldBe = f->Be;
  f->cnBe = f->Be;
  f->oldPs = f->Ps;
  f->cnPs = f->Ps;
  f->oldXi = f->Xi;
  f->cnXi = f->Xi;
  f->oldPi = f->Pi;
  f->cnPi = f->Pi;

  
  return;
}

int solve_t0_slow(const VD& f_xi, const VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
		  PAR *p, int lastpt)
{
  lapack_int N_0 = 3*(lastpt + 1);
  lapack_int kl = 2;
  lapack_int ku = 2;
  lapack_int nrhs = 1;
  lapack_int ldab = 2*kl + ku + 1;
  lapack_int ldb_0 = N_0;
  vector<lapack_int> ipiv_0(N_0);
  VD res_0(ldb_0, 0);
  lapack_int info = 0;
  dbl ell_tol = 20 * p->ell_tol;
  int ell_maxit = p->maxit * 100;

  get_ell_res_abpclean_join(res_0, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt);
  dbl res = norm_inf(res_0);
  int ell_itn = 0;
  while (res > ell_tol) {
    VD jac(ldab*N_0, 0);
    set_jacCMabpslow(jac, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt + 1, kl, ku, ldab);
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N_0, kl, ku, nrhs,
			 &jac[0], ldab, &ipiv_0[0], &res_0[0], ldb_0);
    if (info != 0) { cout << ell_itn << "\nERROR: cannot solve initial elliptics\ninfo = " << info << endl; }
    
    apply_up_abp_join(res_0, f_al, f_be, f_ps, p, lastpt + 1);
    get_ell_res_abpclean_join(res_0, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt);
    res = norm_inf(res_0);

    if (outgoing_null_b(f_al, f_be, f_ps, p, lastpt) <= 0) {
      cout << "\nt = 0\nitn = " << ell_itn << "\nres = " << res << endl;
      return -(lastpt);
    }
    int k = lastpt;
    while (--k > 0) {
      if (outgoing_null(f_al, f_be, f_ps, p, k) <= 0) {
	cout << "\nt = 0\nitn = " << ell_itn << "\nres = " << res << endl;
	return -k;
      }
    }
    if (outgoing_null_f(f_al, f_be, f_ps, p, 0) <= 0) {
      cout << "\nt = 0\nitn = " << ell_itn << "\nres = " << res << endl;
      return -(lastpt + 1);
    }
    if (++ell_itn > ell_maxit) { 
      cout << "\nres = " << res << endl;
      if (res < p->tol) {
	return ell_maxit;
      }
      else { return -2*lastpt; }
    }
    
  }
  return ell_itn;
}



int solve_E_fast(VD& old_xi, VD& old_pi, VD& old_al, VD& old_be, VD& old_ps,
		 VD& f_xi, VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
		 VD& cn_xi, VD& cn_pi, VD& cn_al, VD& cn_be, VD& cn_ps,
		 VD& res_hyp, VD& res_ell, const VD& jac_zero, PAR *p, MAPID& r,
		 int lastpt, int maxit, int i,
		 int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv, int ldb)
{
  string error_response = "idc";
  // set old_f = f = cn_f
  old_xi = f_xi; old_pi = f_pi; old_al = f_al; old_be = f_be; old_ps = f_ps;
  cn_xi = f_xi; cn_pi = f_pi; cn_al = f_al; cn_be = f_be; cn_ps = f_ps;

  VD jac = jac_zero;
  lapack_int info = 0;
  int itn = 0, hyp_itn, ell_itn;
  dbl res = (p->tol) + 1;
  while (res > (p->tol)) {
    hyp_itn = 0; ell_itn = 0;
    while (res > (p->tol)) {
      hyp_solve_px_fast(old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps,
			p, p->r, lastpt);
      res = get_hyp_res_fast(res_hyp, old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps,
			     p, p->r, lastpt);
      if (++hyp_itn > maxit) {
	cout << endl << i << " hyperbolic solver STUCK at t = " << i*(p->dt) << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }
    
    get_ell_res_abp_fast(res_ell, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt);
    res = norm_inf(res_ell);
    while (res > (p->ell_tol)) {
      jac = jac_zero;
      set_jacCMabpfast(jac, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt + 1, kl, ku, ldab);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, nrhs,
			   &jac[0], ldab, &ipiv[0], &res_ell[0], ldb);
      if (info != 0) { cout << i*(p->dt) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl; }
      apply_up_abp_join(res_ell, f_al, f_be, f_ps, p, lastpt + 1);
      get_ell_res_abp_fast(res_ell, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt);
      res = norm_inf(res_ell);
      if (++ell_itn > maxit) {
	cout << endl << i << " elliptic solver STUCK at t = " << i*(p->dt) << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }
    set3_cn(old_ps, old_be, old_al, f_ps, f_be, f_al, cn_ps, cn_be, cn_al, lastpt + 1);
    res = get_hyp_res_fast(res_hyp, old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps,
			   p, p->r, lastpt);
    if (++itn > maxit) {
      cout << endl << i << " solver STUCK at t = " << i*(p->dt) << "\nres = " << res << endl;
      cout << endl << "continue? " << endl;
      cin >> error_response;
    }   
  }  
  return itn;
}

int solve_E_slow(VD& old_xi, VD& old_pi, VD& old_al, VD& old_be, VD& old_ps,
		 VD& f_xi, VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
		 VD& cn_xi, VD& cn_pi, VD& cn_al, VD& cn_be, VD& cn_ps,
		 VD& res_hyp, VD& res_ell, const VD& jac_zero, PAR *p, MAPID& r,
		 int lastpt, int maxit, int i,
		 int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv, int ldb)
{
  string error_response = "idc";
  // set old_f = f = cn_f
  old_xi = f_xi; old_pi = f_pi; old_al = f_al; old_be = f_be; old_ps = f_ps;
  cn_xi = f_xi; cn_pi = f_pi; cn_al = f_al; cn_be = f_be; cn_ps = f_ps;

  VD jac = jac_zero;
  lapack_int info = 0;
  int itn = 0, hyp_itn, ell_itn;
  dbl res = (p->tol) + 1;
  while (res > (p->tol)) {
    hyp_itn = 0; ell_itn = 0;
    while (res > (p->tol)) {
      hyp_solve_px(old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps,
		   p, p->r, lastpt);
      res = get_hyp_res(res_hyp, old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps,
			p, p->r, lastpt);
      if (++hyp_itn > maxit) {
	cout << endl << i << " hyperbolic solver STUCK at t = " << i*(p->dt) << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }
    
    get_ell_res_abpclean_join(res_ell, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt);
    res = norm_inf(res_ell);
    while (res > (p->ell_tol)) {
      jac = jac_zero;
      set_jacCMabpslow(jac, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt + 1, kl, ku, ldab);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, nrhs,
			   &jac[0], ldab, &ipiv[0], &res_ell[0], ldb);
      if (info != 0) { cout << i*(p->dt) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl; }
      apply_up_abp_join(res_ell, f_al, f_be, f_ps, p, lastpt + 1);
      get_ell_res_abpclean_join(res_ell, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt);
      res = norm_inf(res_ell);
      if (++ell_itn > maxit) {
	cout << endl << i << " elliptic solver STUCK at t = " << i*(p->dt) << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }
    set3_cn(old_ps, old_be, old_al, f_ps, f_be, f_al, cn_ps, cn_be, cn_al, lastpt + 1);
    res = get_hyp_res(res_hyp, old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps,
		      p, p->r, lastpt);
    if (++itn > maxit) {
      cout << endl << i << " solver STUCK at t = " << i*(p->dt) << "\nres = " << res << endl;
      cout << endl << "continue? " << endl;
      cin >> error_response;
    }   
  }  
  return itn;
}

int solve_H_slow(VD& old_xi, VD& old_pi, VD& old_al, VD& old_be, VD& old_ps,
		 VD& f_xi, VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
		 VD& cn_xi, VD& cn_pi, VD& cn_al, VD& cn_be, VD& cn_ps,
		 VD& res_hyp, VD& res_ell, const VD& jac_zero, PAR *p, MAPID& r,
		 int lastpt, int maxit, int i,
		 int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv, int ldb)
{
  string error_response = "idc";
  // set old_f = f = cn_f
  old_xi = f_xi; old_pi = f_pi; old_al = f_al; old_be = f_be; old_ps = f_ps;
  cn_xi = f_xi; cn_pi = f_pi; cn_al = f_al; cn_be = f_be; cn_ps = f_ps;

  VD jac = jac_zero;
  lapack_int info = 0;
  int itn = 0, hyp_itn, ell_itn;
  dbl res = (p->tol) + 1;
  while (res > (p->tol)) {
    hyp_itn = 0; ell_itn = 0;
    while (res > (p->tol)) {
      hyp_solve_ppx_slow(old_xi, old_pi, old_ps, f_xi, f_pi, f_ps,
			 cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, p->r, lastpt);
      res = get_hyp_res_ppx(res_hyp, old_xi, old_pi, old_ps, f_xi, f_pi, f_ps,
			    cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, p->r, lastpt);
      if (++hyp_itn > maxit) {
	cout << endl << i << " hyperbolic solver STUCK at t = " << i*(p->dt) << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }
    
    get_ell_res_ab_slow(res_ell, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt);
    res = norm_inf(res_ell);
    while (res > (p->ell_tol)) {
      jac = jac_zero;
      set_jacCM_ab_slow(jac, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt + 1, kl, ku, ldab);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, nrhs,
			   &jac[0], ldab, &ipiv[0], &res_ell[0], ldb);
      if (info != 0) { cout << i*(p->dt) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl; }
      apply_up_ab(res_ell, f_al, f_be, p, lastpt + 1);
      get_ell_res_ab_slow(res_ell, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt);
      res = norm_inf(res_ell);
      if (++ell_itn > maxit) {
	cout << endl << i << " elliptic solver STUCK at t = " << i*(p->dt) << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }
    set2_cn(old_be, old_al, f_be, f_al, cn_be, cn_al, lastpt + 1);
    res = get_hyp_res_ppx(res_hyp, old_xi, old_pi, old_ps, f_xi, f_pi, f_ps,
			  cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, p->r, lastpt);
    if (++itn > maxit) {
      cout << endl << i << " solver STUCK at t = " << i*(p->dt) << "\nres = " << res << endl;
      cout << endl << "continue? " << endl;
      cin >> error_response;
    }   
  }  
  return itn;
}

int solve_Hptol_slow(VD& old_xi, VD& old_pi, VD& old_al, VD& old_be, VD& old_ps,
		     VD& f_xi, VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
		     VD& cn_xi, VD& cn_pi, VD& cn_al, VD& cn_be, VD& cn_ps,
		     VD& res_hyp, VD& res_ell, const VD& jac_zero, PAR *p, MAPID& r,
		     int lastpt, int maxit, int i,
		     int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv, int ldb)
{
  string error_response = "idc";
  // set old_f = f = cn_f
  old_xi = f_xi; old_pi = f_pi; old_al = f_al; old_be = f_be; old_ps = f_ps;
  cn_xi = f_xi; cn_pi = f_pi; cn_al = f_al; cn_be = f_be; cn_ps = f_ps;

  VD jac = jac_zero;
  lapack_int info = 0;
  int itn = 0, hyp_itn, ell_itn;
  dbl res = (p->tol) + 1;

  int kps = 2 * (lastpt + 1);
  while (res > (p->tol)) {
    hyp_itn = 0; ell_itn = 0;
    while (res > (p->ell_tol)) {
      hyp_solve_PSI_ONLY(old_ps, f_ps, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, p->r, lastpt);      
      res_hyp[kps] = neumann0res(f_ps, p);
      for (int k = 1; k < lastpt; ++k) {
	res_hyp[kps + k] = fda_hyp_resPs(old_ps, f_ps, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
      }
      res_hyp[kps + lastpt] = fdaR_hyp_resPs(f_ps, p, lastpt);
      res =  max(  *max_element(res_hyp.begin() + kps, res_hyp.end()),
		   -(*min_element(res_hyp.begin() + kps, res_hyp.end()))  );
      if (++hyp_itn > maxit) {
	cout << endl << i << " PSI hyperbolic solver STUCK at t = " << i*(p->dt) << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }
    res = (p->tol) + 1;
    hyp_itn = 0;
    while (res > (p->tol)) {
      hyp_solve_px(old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, p->r, lastpt);
      res = get_hyp_res(res_hyp, old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, p->r, lastpt);
      if (++hyp_itn > maxit) {
	cout << endl << i << " hyperbolic solver STUCK at t = " << i*(p->dt) << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }    
    get_ell_res_ab_slow(res_ell, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt);
    res = norm_inf(res_ell);
    while (res > (p->ell_tol)) {
      jac = jac_zero;
      set_jacCM_ab_slow(jac, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt + 1, kl, ku, ldab);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, nrhs,
			   &jac[0], ldab, &ipiv[0], &res_ell[0], ldb);
      if (info != 0) { cout << i*(p->dt) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl; }
      apply_up_ab(res_ell, f_al, f_be, p, lastpt + 1);
      get_ell_res_ab_slow(res_ell, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt);
      res = norm_inf(res_ell);
      if (++ell_itn > maxit) {
	cout << endl << i << " elliptic solver STUCK at t = " << i*(p->dt) << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }
    set2_cn(old_al, old_be, f_al, f_be, cn_al, cn_be, lastpt + 1);
    for (int k = 1; k < lastpt; ++k) {
      res_hyp[kps + k] = fda_hyp_resPs(old_ps, f_ps, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    }
    res =  max(  *max_element(res_hyp.begin() + kps, res_hyp.end()),
		   -(*min_element(res_hyp.begin() + kps, res_hyp.end()))  );
    if (res < (p->ell_tol)) { 
      res = get_hyp_res(res_hyp, old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, p->r, lastpt);
    }
    else { res = (p->tol) + 1; }
    if (++itn > maxit) {
      cout << endl << i << " solver STUCK at t = " << i*(p->dt) << "\nres = " << res << endl;
      cout << endl << "continue? " << endl;
      cin >> error_response;
    }   
  }  
  return itn;
}

int solve_Hsearch(VD& old_xi, VD& old_pi, VD& old_al, VD& old_be, VD& old_ps,
		  VD& f_xi, VD& f_pi, VD& f_al, VD& f_be, VD& f_ps,
		  VD& cn_xi, VD& cn_pi, VD& cn_al, VD& cn_be, VD& cn_ps,
		  VD& res_hyp, VD& res_ell, const VD& jac_zero, PAR *p, MAPID& r,
		  int lastpt, int maxit, int i,
		  int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv, int ldb)
{
  string error_response = "idc";
  // set old_f = f = cn_f
  old_xi = f_xi; old_pi = f_pi; old_al = f_al; old_be = f_be; old_ps = f_ps;
  cn_xi = f_xi; cn_pi = f_pi; cn_al = f_al; cn_be = f_be; cn_ps = f_ps;

  VD jac = jac_zero;
  lapack_int info = 0;
  int itn = 0, hyp_itn, ell_itn;
  dbl res = (p->tol) + 1;

  int kps = 2 * (lastpt + 1);
  while (res > (p->tol)) {
    hyp_itn = 0; ell_itn = 0;
    while (res > (p->ell_tol)) {
      hyp_solve_PSI_ONLY(old_ps, f_ps, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, p->r, lastpt);      
      res_hyp[kps] = neumann0res(f_ps, p);
      for (int k = 1; k < lastpt; ++k) {
	res_hyp[kps + k] = fda_hyp_resPs(old_ps, f_ps, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
      }
      res_hyp[kps + lastpt] = fdaR_hyp_resPs(f_ps, p, lastpt);
      res =  max(  *max_element(res_hyp.begin() + kps, res_hyp.end()),
		   -(*min_element(res_hyp.begin() + kps, res_hyp.end()))  );
      if (++hyp_itn > maxit) {
	cout << endl << i << " PSI hyperbolic solver STUCK at t = " << i*(p->dt) << endl;
	cout << endl << "continue? " << endl;
	cin >> error_response;
      }
    }
    res = (p->tol) + 1;
    hyp_itn = 0;
    while (res > (p->tol)) {
      hyp_solve_px(old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, p->r, lastpt);
      res = get_hyp_res(res_hyp, old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, p->r, lastpt);
      if (++hyp_itn > maxit) {
	cout << endl << i << " hyperbolic solver STUCK at t = " << i*(p->dt) << endl;
        cout << "res = " << res << endl;
        return -1;
      }
    }    
    get_ell_res_ab_slow(res_ell, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt);
    res = norm_inf(res_ell);
    while (res > (p->ell_tol)) {
      jac = jac_zero;
      set_jacCM_ab_slow(jac, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt + 1, kl, ku, ldab);
      info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, nrhs,
			   &jac[0], ldab, &ipiv[0], &res_ell[0], ldb);
      if (info != 0) { cout << i*(p->dt) << "\nERROR: cannot solve elliptic equations\ninfo = " << info << endl; }
      apply_up_ab(res_ell, f_al, f_be, p, lastpt + 1);
      get_ell_res_ab_slow(res_ell, f_xi, f_pi, f_al, f_be, f_ps, p, p->r, lastpt);
      res = norm_inf(res_ell);
      if (++ell_itn > maxit) {
	cout << endl << i << " elliptic solver STUCK at t = " << i*(p->dt) << endl;
	cout << "res = " << res << endl;
	if (res < (p->tol)) { res = 0; }
        else { return -2; }
      }
    }
    set2_cn(old_al, old_be, f_al, f_be, cn_al, cn_be, lastpt + 1);
    for (int k = 1; k < lastpt; ++k) {
      res_hyp[kps + k] = fda_hyp_resPs(old_ps, f_ps, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    }
    res =  max(  *max_element(res_hyp.begin() + kps, res_hyp.end()),
		   -(*min_element(res_hyp.begin() + kps, res_hyp.end()))  );
    if (res < (p->ell_tol)) { 
      res = get_hyp_res(res_hyp, old_xi, old_pi, f_xi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, p->r, lastpt);
    }
    else { res = (p->tol) + 1; }
    if (++itn > maxit) {
      cout << endl << i << " solver STUCK at t = " << i*(p->dt) << "\nres = " << res << endl;
      cout << "res = " << res << endl;
      return -3;
    }   
  }  
  return itn;
}

#endif
