#ifndef EKG_CLEAN_H_INCLUDED
#define EKG_CLEAN_H_INCLUDED

#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <vector> // for everything
#include <cmath> // for ICs
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf
#include "sim-header.h"
#include "sim-structs.h"
#include "fda-io.h"
#include "fda-fns.h"
#include "ekg-fns.h"
#include "jacobian.h"
#include "ekg-proc.h"
#include "lapacke.h"

using namespace std;

// perform gauss-seidel update on xi, pi (xpp version for including f_ps)
void hyp_solve_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
		  const VD& cn_al, const VD& cn_be, const VD& cn_ps, PAR *p, MAPID& r, int lastpt);
void hyp_solve_px(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
		  const VD& cn_al, const VD& cn_be, const VD& cn_ps, PAR *p, MAPID& r, int lastpt);
void hyp_solve_px_fast(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
		       const VD& cn_al, const VD& cn_be, const VD& cn_ps, PAR *p, MAPID& r, int lastpt);
void hyp_solve(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
	       const VD& cn_al, const VD& cn_be, const VD& cn_ps, PAR *p, MAPID& r, int lastpt);
void hyp_solve_ppx_slow(const VD& old_xi, const VD& old_pi, const VD& old_ps,
			VD& f_xi, VD& f_pi, VD& f_ps, VD& cn_xi, VD& cn_pi,
			const VD& cn_al, const VD& cn_be, const VD& cn_ps, PAR *p, MAPID& r, int lastpt);
void hyp_solve_PSI_ONLY(const VD& old_ps, VD& f_ps, const VD& cn_xi, const VD& cn_pi, const VD& cn_al,
			const VD& cn_be, VD& cn_ps, PAR *p, MAPID& r, int lastpt);
// set res_vals to be residual = L[f] for rhs of jacobian.(-delta) = residual
// get inf-norm of residuals, residuals vector set
dbl full_res0clean(VD& residuals, const VD& old_xi, const VD& old_pi,
		   const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		   const VD& cn_xi, const VD& cn_pi, const VD& cn_al, const VD& cn_be, const VD& cn_ps,
		   PAR *p, MAPID& r, int lastpt);
dbl get_hyp_res(VD& residuals, const VD& old_xi, const VD& old_pi, const VD& f_xi, const VD& f_pi,
		const VD& cn_xi, const VD& cn_pi, const VD& cn_al, const VD& cn_be, const VD& cn_ps,
		PAR *p, MAPID& r, int lastpt);
dbl get_hyp_res_ppx(VD& res_hyp, const VD& old_xi, const VD& old_pi, const VD& old_ps,
		    const VD& f_xi, const VD& f_pi, const VD& f_ps,
		    const VD& cn_xi, const VD& cn_pi, const VD& cn_al, const VD& cn_be, const VD& cn_ps,
		    PAR *p, MAPID& r, int lastpt);
dbl get_hyp_res_fast(VD& residuals, const VD& old_xi, const VD& old_pi, const VD& f_xi, const VD& f_pi,
		     const VD& cn_xi, const VD& cn_pi, const VD& cn_al, const VD& cn_be, const VD& cn_ps,
		     PAR *p, MAPID& r, int lastpt);
void get_ell_res_abpclean_join(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
			       const VD& f_ps, PAR *p, MAPID& r, int lastpt);
void get_ell_res_abp_fast(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
			  const VD& f_ps, PAR *p, MAPID& r, int lastpt);
void get_ell_res_ab_slow(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
			 const VD& f_ps, PAR *p, MAPID& r, int lastpt);
void get_ell_res_pbaclean_join(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
			       const VD& f_ps, PAR *p, MAPID& r, int lastpt);
double get_ell_res_ind(VD& res_al, VD& res_be, VD& res_ps, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
		       const VD& f_ps, PAR *p, MAPID& r, int lastpt);
inline void get_resAl(VD& res_al, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		      PAR *p, MAPID& r, int lastpt);
inline void get_resBe(VD& res_be, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		      PAR *p, MAPID& r, int lastpt);
inline void get_resPs(VD& res_ps, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		      PAR *p, MAPID& r, int lastpt);
void get_resPs_hyp(VD& res_ps, const VD& old_ps, const VD& f_ps, const VD& cn_xi, const VD& cn_pi,
		   const VD& cn_al, const VD& cn_be, const VD& cn_ps, PAR *p, MAPID& r, int lastwr, int save_pt);
void get_resPs_ell(VD& res_ps, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		   PAR *p, MAPID& r, int lastwr, int save_pt);
inline void apply_up_field(const VD& deltas, VD& field, PAR *p, int npts);
inline void apply_up_ab(const VD& deltas, VD& f_al, VD& f_be, PAR *p, int npts);
inline void apply_up_ab_cn(const VD& deltas, VD& f_al, VD& f_be, PAR *p, int npts);
inline void apply_up_abp_join(const VD& deltas, VD& f_al, VD& f_be, VD& f_ps, PAR *p, int npts);
inline void apply_up_abp_join_cn(const VD& deltas, const VD& old_al, const VD& old_be, const VD& old_ps,
				 VD& f_al, VD& f_be, VD& f_ps, VD& cn_al, VD& cn_be, VD& cn_ps, PAR *p, int npts);
int search_for_horizon(const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, MAPID& r, int lastpt, int i, double t)
{
  string horizon_response = "awesome";
  if (outgoing_null_b(f_al, f_be, f_ps, p, lastpt) <= 0) {
    cout << i << " at t = " << t << "\nAPPARENT HORIZON FOUND at r = " << (p->rmax) << "\nk = " << lastpt;
    cout << "\nr_areal = " << (p->rmax)*sq(f_ps[lastpt]) << endl;
    return 1;
  }
  int k = lastpt;
  while (--k > 0) {
    if (outgoing_null(f_al, f_be, f_ps, p, k) <= 0) {
      cout << i << " at t = " << t << "\nAPPARENT HORIZON FOUND at r = " << r[k] << "\nk = " << k;
      cout << "\nr_areal = " << r[k]*sq(f_ps[k]) << endl;
      return 1;
    }
  }
  if (outgoing_null_f(f_al, f_be, f_ps, p, 0) <= 0) {
    cout << i << " at t = " << t << "\nAPPARENT HORIZON FOUND at r = " << (p->rmin) << "\nk = " << 0;
    cout << "\nr_areal = " << (p->rmin)*sq(f_ps[0]) << endl;
    return 1;
  }
  return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// ***************************** HYP_SOLVE
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// perform gauss-seidel update (xi and pi only)
////////////////////////////////////////////////////////////////////////////////////////////////
void hyp_solve_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
		  const VD& cn_al, const VD& cn_be, const VD& cn_ps, PAR *p, MAPID& r, int lastpt)
{
  neumann0(f_pi);
  cn_pi[0] = 0.5 * (old_pi[0] + f_pi[0]);
  for (int k = 1; k < lastpt; ++k) {
    // update f_xi & cn_xi
    f_xi[k] = fda_xi(old_xi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    cn_xi[k] = 0.5 * (old_xi[k] + f_xi[k]);
    // update f_pi & cn_pi
    f_pi[k] = fda_pi(old_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    cn_pi[k] = 0.5 * (old_pi[k] + f_pi[k]);
  }
  // r = R BOUNDARY
  sommerfeld(old_xi, f_xi, p, lastpt);
  cn_xi[lastpt] = 0.5 * (old_xi[lastpt] + f_xi[lastpt]);
  sommerfeld(old_pi, f_pi, p, lastpt);
  cn_pi[lastpt] = 0.5 * (old_pi[lastpt] + f_pi[lastpt]);
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
void hyp_solve_px(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
		  const VD& cn_al, const VD& cn_be, const VD& cn_ps, PAR *p, MAPID& r, int lastpt)
{
  neumann0(f_pi);
  cn_pi[0] = 0.5 * (old_pi[0] + f_pi[0]);
  for (int k = 1; k < lastpt; ++k) {
    // update f_pi & cn_pi
    f_pi[k] = fda_pi(old_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    cn_pi[k] = 0.5 * (old_pi[k] + f_pi[k]);
    // update f_xi & cn_xi
    f_xi[k] = fda_xi(old_xi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    cn_xi[k] = 0.5 * (old_xi[k] + f_xi[k]);
  }
  // r = R BOUNDARY
  sommerfeld(old_pi, f_pi, p, lastpt);
  cn_pi[lastpt] = 0.5 * (old_pi[lastpt] + f_pi[lastpt]);
  sommerfeld(old_xi, f_xi, p, lastpt);
  cn_xi[lastpt] = 0.5 * (old_xi[lastpt] + f_xi[lastpt]);
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
void hyp_solve_px_fast(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
		       const VD& cn_al, const VD& cn_be, const VD& cn_ps, PAR *p, MAPID& r, int lastpt)
{
  double lam6_part;
  double ps2_m, ps2 = sq(cn_ps[0]), ps2_p = sq(cn_ps[1]);
  double al_ps2_m, al_ps2 = cn_al[0] / ps2, al_ps2_p = cn_al[1] / ps2_p;
  neumann0(f_pi);
  cn_pi[0] = 0.5 * (old_pi[0] + f_pi[0]);
  for (int k = 1; k < lastpt; ++k) {
    lam6_part = r[LAM6VAL] * (d_c(cn_be,k) + cn_be[k]*(4*r[DRVAL]*r[-k] + 6*d_c(cn_ps,k)/cn_ps[k]));
    ps2_m = ps2; ps2 = ps2_p; ps2_p = sq(cn_ps[k+1]);
    al_ps2_m = al_ps2; al_ps2 = al_ps2_p; al_ps2_p = cn_al[k+1] / ps2_p;
    // update f_pi & cn_pi
    f_pi[k] = ( old_pi[k]*(1 - lam6_part) + (r[LAM2VAL] / sq(r[k]*ps2)) *
		(sq(r[k+1]*ps2_p)*(al_ps2_p*cn_xi[k+1] + cn_be[k+1]*cn_pi[k+1])
		 - sq(r[k-1]*ps2_m)*(al_ps2_m*cn_xi[k-1] + cn_be[k-1]*cn_pi[k-1])) )
      / (1 + lam6_part);
    cn_pi[k] = 0.5 * (old_pi[k] + f_pi[k]);
    // update f_xi & cn_xi
    f_xi[k] = old_xi[k] + r[LAM2VAL]*(al_ps2_p*cn_pi[k+1] + cn_be[k+1]*cn_xi[k+1]
				      - al_ps2_m*cn_pi[k-1] - cn_be[k-1]*cn_xi[k-1]);
    cn_xi[k] = 0.5 * (old_xi[k] + f_xi[k]);
  }
  // r = R BOUNDARY
  sommerfeld(old_pi, f_pi, p, lastpt);
  cn_pi[lastpt] = 0.5 * (old_pi[lastpt] + f_pi[lastpt]);
  sommerfeld(old_xi, f_xi, p, lastpt);
  cn_xi[lastpt] = 0.5 * (old_xi[lastpt] + f_xi[lastpt]);
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
void hyp_solve(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi, VD& cn_xi, VD& cn_pi,
	       const VD& cn_al, const VD& cn_be, const VD& cn_ps, PAR *p, MAPID& r, int lastpt)
{
  neumann0(f_pi);
  cn_pi[0] = 0.5 * (old_pi[0] + f_pi[0]);
  for (int k = 1; k < lastpt; ++k) {
    // update f_pi & f_xi
    f_pi[k] = fda_pi(old_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    f_xi[k] = fda_xi(old_xi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    // update cn_pi & cn_xi
    cn_pi[k] = 0.5 * (old_pi[k] + f_pi[k]);
    cn_xi[k] = 0.5 * (old_xi[k] + f_xi[k]);
  }
  // r = R BOUNDARY
  sommerfeld(old_pi, f_pi, p, lastpt);
  sommerfeld(old_xi, f_xi, p, lastpt);
  cn_pi[lastpt] = 0.5 * (old_pi[lastpt] + f_pi[lastpt]);
  cn_xi[lastpt] = 0.5 * (old_xi[lastpt] + f_xi[lastpt]);
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
void hyp_solve_ppx_slow(const VD& old_xi, const VD& old_pi, const VD& old_ps, VD& f_xi, VD& f_pi, VD& f_ps,
			VD& cn_xi, VD& cn_pi, const VD& cn_al, const VD& cn_be, VD& cn_ps, PAR *p, MAPID& r, int lastpt)
{
  neumann0(f_ps);
  cn_ps[0] = 0.5 * (old_ps[0] + f_ps[0]);
  neumann0(f_pi);
  cn_pi[0] = 0.5 * (old_pi[0] + f_pi[0]);
  for (int k = 1; k < lastpt; ++k) {
    // update f_ps & cn_ps
    f_ps[k] = fda_hyp_ps(old_ps, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    cn_ps[k] = 0.5 * (old_ps[k] + f_ps[k]);
    // update f_pi & cn_pi
    f_pi[k] = fda_pi(old_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    cn_pi[k] = 0.5 * (old_pi[k] + f_pi[k]);
    // update f_xi & cn_xi
    f_xi[k] = fda_xi(old_xi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    cn_xi[k] = 0.5 * (old_xi[k] + f_xi[k]);
  }
  // r = R BOUNDARY
  f_ps[lastpt] = fdaR_hyp_ps(f_ps, p, lastpt);
  cn_ps[lastpt] = 0.5 * (old_ps[lastpt] + f_ps[lastpt]);
  sommerfeld(old_pi, f_pi, p, lastpt);
  cn_pi[lastpt] = 0.5 * (old_pi[lastpt] + f_pi[lastpt]);
  sommerfeld(old_xi, f_xi, p, lastpt);
  cn_xi[lastpt] = 0.5 * (old_xi[lastpt] + f_xi[lastpt]);
  return;
}

void hyp_solve_PSI_ONLY(const VD& old_ps, VD& f_ps, const VD& cn_xi, const VD& cn_pi, const VD& cn_al,
			const VD& cn_be, VD& cn_ps, PAR *p, MAPID& r, int lastpt)
{
  neumann0(f_ps);
  cn_ps[0] = 0.5 * (old_ps[0] + f_ps[0]);
  for (int k = 1; k < lastpt; ++k) {
    // update f_ps & cn_ps
    f_ps[k] = fda_hyp_ps(old_ps, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    cn_ps[k] = 0.5 * (old_ps[k] + f_ps[k]);
  }
  // r = R BOUNDARY
  f_ps[lastpt] = fdaR_hyp_ps(f_ps, p, lastpt);
  cn_ps[lastpt] = 0.5 * (old_ps[lastpt] + f_ps[lastpt]);
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// set res_vals to be residual = L[f] for rhs of jacobian.(-delta) = residual
// ************************  GET_FULL_RES
////////////////////////////////////////////////////////////////////////////////////////////////
dbl full_res0clean(VD& residuals, const VD& old_xi, const VD& old_pi,
		   const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		   const VD& cn_xi, const VD& cn_pi, const VD& cn_al, const VD& cn_be, const VD& cn_ps,
		   PAR *p, MAPID& r, int lastpt)
{
  int kpi = lastpt + 1, kal = 2*kpi, kbe = 3*kpi, kps = 4*kpi;
  residuals[0] = dirichlet0res(f_xi);
  residuals[kpi] = neumann0res(f_pi, p);
  residuals[kal] = neumann0res(f_al, p);
  residuals[kbe] = dirichlet0res(f_be);
  residuals[kps] = neumann0res(f_ps, p);
  for (int k = 1; k < lastpt; ++k) {
    residuals[k] = fda_resXi(old_xi, f_xi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    residuals[kpi+k] = fda_resPi(old_pi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    residuals[kal+k] = fda_resAl(f_xi, f_pi, f_al, f_be, f_ps, p, k);
    residuals[kbe+k] = fda_resBe(f_xi, f_pi, f_al, f_be, f_ps, p, k);
    residuals[kps+k] = fda_resPs(f_xi, f_pi, f_al, f_be, f_ps, p, k);
  }
  residuals[lastpt] = sommerfeldres(old_xi, f_xi, p, lastpt);
  residuals[kpi] = sommerfeldres(old_pi, f_pi, p, lastpt);
  residuals[kal] = fdaR_resAl(f_al, p, lastpt);
  residuals[kbe] = fdaR_resBe(f_be, p, lastpt);
  residuals[kps] = fdaR_resPs(f_ps, p, lastpt);
  return max(  *max_element(residuals.begin(), residuals.end()),
	     -(*min_element(residuals.begin(), residuals.end()))  );
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// ***************************** GET_HYP_RES
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

dbl get_hyp_res(VD& res_hyp, const VD& old_xi, const VD& old_pi, const VD& f_xi, const VD& f_pi,
		const VD& cn_xi, const VD& cn_pi, const VD& cn_al, const VD& cn_be, const VD& cn_ps,
		PAR *p, MAPID& r, int lastpt)
{
  int kpi = lastpt + 1;
  res_hyp[0] = dirichlet0res(f_xi);
  res_hyp[kpi] = neumann0res(f_pi, p);
  for (int k = 1; k < lastpt; ++k) {
    res_hyp[k] = fda_resXi(old_xi, f_xi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    res_hyp[kpi + k] = fda_resPi(old_pi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
  }
  res_hyp[lastpt] = sommerfeldres(old_xi, f_xi, p, lastpt);
  res_hyp[kpi + lastpt] = sommerfeldres(old_pi, f_pi, p, lastpt);
  return max(  *max_element(res_hyp.begin(), res_hyp.end()),
	      -(*min_element(res_hyp.begin(), res_hyp.end()))  );
}

dbl get_hyp_res_ppx(VD& res_hyp, const VD& old_xi, const VD& old_pi, const VD& old_ps,
		    const VD& f_xi, const VD& f_pi, const VD& f_ps,
		    const VD& cn_xi, const VD& cn_pi, const VD& cn_al, const VD& cn_be, const VD& cn_ps,
		    PAR *p, MAPID& r, int lastpt)
{
  int kpi = lastpt + 1;
  int kps = 2 * kpi;
  res_hyp[0] = dirichlet0res(f_xi);
  res_hyp[kpi] = neumann0res(f_pi, p);
  res_hyp[kps] = neumann0res(f_ps, p);
  for (int k = 1; k < lastpt; ++k) {
    res_hyp[k] = fda_resXi(old_xi, f_xi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    res_hyp[kpi + k] = fda_resPi(old_pi, f_pi, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
    res_hyp[kps + k] = fda_hyp_resPs(old_ps, f_ps, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, k);
  }
  res_hyp[lastpt] = sommerfeldres(old_xi, f_xi, p, lastpt);
  res_hyp[kpi + lastpt] = sommerfeldres(old_pi, f_pi, p, lastpt);
  res_hyp[kps + lastpt] = fdaR_hyp_resPs(f_ps, p, lastpt);
  return max(  *max_element(res_hyp.begin(), res_hyp.end()),
	      -(*min_element(res_hyp.begin(), res_hyp.end()))  );
}

dbl get_hyp_res_fast(VD& res_hyp, const VD& old_xi, const VD& old_pi, const VD& f_xi, const VD& f_pi,
		     const VD& cn_xi, const VD& cn_pi, const VD& cn_al, const VD& cn_be, const VD& cn_ps,
		     PAR *p, MAPID& r, int lastpt)
{
  double lam6_part;
  double ps2_m, ps2 = sq(cn_ps[0]), ps2_p = sq(cn_ps[1]);
  double al_ps2_m, al_ps2 = cn_al[0] / ps2, al_ps2_p = cn_al[1] / ps2_p;
  int kpi = lastpt + 1;
  res_hyp[0] = dirichlet0res(f_xi);
  res_hyp[kpi] = neumann0res(f_pi, p);
  for (int k = 1; k < lastpt; ++k) {
    lam6_part = r[LAM6VAL] * (d_c(cn_be,k) + cn_be[k]*(4*r[DRVAL]*r[-k] + 6*d_c(cn_ps,k)/cn_ps[k]));
    ps2_m = ps2; ps2 = ps2_p; ps2_p = sq(cn_ps[k+1]);
    al_ps2_m = al_ps2; al_ps2 = al_ps2_p; al_ps2_p = cn_al[k+1] / ps2_p;
    
    res_hyp[k] = r[INDT]*( f_xi[k] - old_xi[k] -
			   r[LAM2VAL]*(al_ps2_p*cn_pi[k+1] + cn_be[k+1]*cn_xi[k+1]
				       - al_ps2_m*cn_pi[k-1] - cn_be[k-1]*cn_xi[k-1]) );
    res_hyp[kpi + k] = r[INDT]*( f_pi[k]*(lam6_part + 1) + old_pi[k]*(lam6_part - 1)
				 - (r[LAM2VAL] / sq(r[k]*ps2)) *
				 (sq(r[k+1]*ps2_p)*(al_ps2_p*cn_xi[k+1] + cn_be[k+1]*cn_pi[k+1])
				  - sq(r[k-1]*ps2_m)*(al_ps2_m*cn_xi[k-1] + cn_be[k-1]*cn_pi[k-1])) );
  }
  res_hyp[lastpt] = sommerfeldres(old_xi, f_xi, p, lastpt);
  res_hyp[kpi + lastpt] = sommerfeldres(old_pi, f_pi, p, lastpt);
  return max(  *max_element(res_hyp.begin(), res_hyp.end()),
	      -(*min_element(res_hyp.begin(), res_hyp.end()))  );
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// ***************************** GET_ELL_RES
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

void get_ell_res_abpclean_join(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
			       const VD& f_ps, PAR *p, MAPID& r, int lastpt)
{
  int kbe = lastpt + 1;
  int kps = 2*kbe;
  res_ell[0] = neumann0res(f_al, p);
  res_ell[kbe] = dirichlet0res(f_be);
  res_ell[kps] = neumann0res(f_ps, p);
  for (int k = 1; k < lastpt; ++k) {
    res_ell[k] = fda_resAl(f_xi, f_pi, f_al, f_be, f_ps, p, k);
    res_ell[kbe + k] = fda_resBe(f_xi, f_pi, f_al, f_be, f_ps, p, k);
    res_ell[kps + k] = fda_resPs(f_xi, f_pi, f_al, f_be, f_ps, p, k);
  }
  res_ell[lastpt] = fdaR_resAl(f_al, p, lastpt);
  res_ell[kbe + lastpt] = fdaR_resBe(f_be, p, lastpt);
  res_ell[kps + lastpt] = fdaR_resPs(f_ps, p, lastpt);
  return;
}

void get_ell_res_abp_fast(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
			  const VD& f_ps, PAR *p, MAPID& r, int lastpt)
{
  double drbe_r, dps_ps, alin, p4db2_a;
  int kbe = lastpt + 1;
  int kps = 2*kbe;
  res_ell[0] = neumann0res(f_al, p);
  res_ell[kbe] = dirichlet0res(f_be);
  res_ell[kps] = neumann0res(f_ps, p);
  for (int k = 1; k < lastpt; ++k) {
    drbe_r = ddr_c(f_be,p,k) - r[-k]*f_be[k];
    dps_ps = ddr_c(f_ps,p,k) / f_ps[k];
    alin = 1 / f_al[k];
    p4db2_a = sq( sq(f_ps[k]) * drbe_r ) * alin; 
    
    
    res_ell[k] = ddr2_c(f_al,p,k) + 2*ddr_c(f_al,p,k)*(r[-k] + dps_ps)
      - r[TWO_THIRDS]*p4db2_a - 8*M_PI*f_al[k]*sq(f_pi[k]);
    
    res_ell[kbe + k] = ddr2_c(f_be,p,k) + drbe_r*(2*r[-k] + 6*dps_ps - ddr_c(f_al,p,k)*alin)
      + 12*M_PI*f_al[k]*f_xi[k]*f_pi[k]/sq(f_ps[k]);
    
    res_ell[kps + k] = ddr2_c(f_ps,p,k) + f_ps[k]*( 2*r[-k]*dps_ps + r[TWELFTH]*p4db2_a*alin
						    + M_PI*(sq(f_xi[k]) + sq(f_pi[k])) );
  }
  res_ell[lastpt] = fdaR_resAl(f_al, p, lastpt);
  res_ell[kbe + lastpt] = fdaR_resBe(f_be, p, lastpt);
  res_ell[kps + lastpt] = fdaR_resPs(f_ps, p, lastpt);
  return;
}

void get_ell_res_ab_slow(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
			 const VD& f_ps, PAR *p, MAPID& r, int lastpt)
{
  int kbe = lastpt + 1;
  res_ell[0] = neumann0res(f_al, p);
  res_ell[kbe] = dirichlet0res(f_be);
  for (int k = 1; k < lastpt; ++k) {
    res_ell[k] = fda_resAl(f_xi, f_pi, f_al, f_be, f_ps, p, k);
    res_ell[kbe + k] = fda_resBe(f_xi, f_pi, f_al, f_be, f_ps, p, k);
  }
  res_ell[lastpt] = fdaR_resAl(f_al, p, lastpt);
  res_ell[kbe + lastpt] = fdaR_resBe(f_be, p, lastpt);
  return;
}

// ******** NOW THIS JOINS P->B->A
void get_ell_res_pbaclean_join(VD& res_ell, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
			       const VD& f_ps, PAR *p, MAPID& r, int lastpt)
{
  int kbe = lastpt + 1;
  int kal = 2*kbe;
  res_ell[0] = neumann0res(f_ps, p);
  res_ell[kbe] = dirichlet0res(f_be);
  res_ell[kal] = neumann0res(f_al, p);
  for (int k = 1; k < lastpt; ++k) {
    res_ell[k] = fda_resPs(f_xi, f_pi, f_al, f_be, f_ps, p, k);
    res_ell[kbe + k] = fda_resBe(f_xi, f_pi, f_al, f_be, f_ps, p, k);
    res_ell[kal + k] = fda_resAl(f_xi, f_pi, f_al, f_be, f_ps, p, k);
  }
  res_ell[lastpt] = fdaR_resPs(f_ps, p, lastpt);
  res_ell[kbe + lastpt] = fdaR_resBe(f_be, p, lastpt);
  res_ell[kal + lastpt] = fdaR_resAl(f_al, p, lastpt);
  return;
}
double get_ell_res_ind(VD& res_al, VD& res_be, VD& res_ps, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be,
		       const VD& f_ps, PAR *p, MAPID& r, int lastpt)
{
  res_ps[0] = neumann0res(f_ps, p);
  res_be[0] = dirichlet0res(f_be);
  res_al[0] = neumann0res(f_al, p);
  for (int k = 1; k < lastpt; ++k) {
    res_ps[k] = fda_resPs(f_xi, f_pi, f_al, f_be, f_ps, p, k);
    res_be[k] = fda_resBe(f_xi, f_pi, f_al, f_be, f_ps, p, k);
    res_al[k] = fda_resAl(f_xi, f_pi, f_al, f_be, f_ps, p, k);
  }
  res_ps[lastpt] = fdaR_resPs(f_ps, p, lastpt);
  res_be[lastpt] = fdaR_resBe(f_be, p, lastpt);
  res_al[lastpt] = fdaR_resAl(f_al, p, lastpt);
  return max(norm_inf(res_ps), max(norm_inf(res_be), norm_inf(res_al)));
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// ***************************** GET_IND_RES
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
inline void get_resAl(VD& res_al, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		      PAR *p, MAPID& r, int lastpt)
{
  res_al[0] = neumann0res(f_al, p);
  for (int k = 0; k < lastpt; ++k) { res_al[k] = fda_resAl(f_xi, f_pi, f_al, f_be, f_ps, p, k); }
  res_al[lastpt] = fdaR_resAl(f_al, p, lastpt);
  return;
}
inline void get_resBe(VD& res_be, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		      PAR *p, MAPID& r, int lastpt)
{
  res_be[0] = dirichlet0res(f_be);
  for (int k = 0; k < lastpt; ++k) { res_be[k] = fda_resBe(f_xi, f_pi, f_al, f_be, f_ps, p, k); }
  res_be[lastpt] = fdaR_resBe(f_al, p, lastpt);
  return;
}
inline void get_resPs(VD& res_ps, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		     PAR *p, MAPID& r, int lastpt)
{
  res_ps[0] = neumann0res(f_ps, p);
  for (int k = 0; k < lastpt; ++k) { res_ps[k] = fda_resPs(f_xi, f_pi, f_al, f_be, f_ps, p, k); }
  res_ps[lastpt] = fdaR_resPs(f_ps, p, lastpt);
  return;
}

void get_resPs_hyp(VD& res_ps, const VD& old_ps, const VD& f_ps, const VD& cn_xi, const VD& cn_pi,
		   const VD& cn_al, const VD& cn_be, const VD& cn_ps, PAR *p, MAPID& r, int lastwr, int save_pt)
{
  res_ps[0] = neumann0res(f_ps, p);
  int s = save_pt;
  for (int k = 0; k < lastwr; ++k) {
    res_ps[k] = fda_hyp_resPs(old_ps, f_ps, cn_xi, cn_pi, cn_al, cn_be, cn_ps, p, s);
    s += save_pt;
  }
  res_ps[lastwr] = fdaR_hyp_resPs(f_ps, p, s);
  return;
}

void get_resPs_ell(VD& res_ps, const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		   PAR *p, MAPID& r, int lastwr, int save_pt)
{
  res_ps[0] = neumann0res(f_ps, p);
  int s = save_pt;
  for (int k = 0; k < lastwr; ++k) {
    res_ps[k] = fda_resPs(f_xi, f_pi, f_al, f_be, f_ps, p, s);
    s += save_pt;
  }
  res_ps[lastwr] = fdaR_resPs(f_ps, p, s);
  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// ***************************** APPLY_UP
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

inline void apply_up_field(const VD& deltas, VD& field, PAR *p, int npts)
{
  for (int k = 0; k < npts; ++k) { field[k] -= (p->ell_up_weight)*deltas[k]; }
  return;
}

inline void apply_up_ab(const VD& deltas, VD& f_al, VD& f_be, PAR *p, int npts)
{
  for (int k = 0; k < npts; ++k) {
    f_al[k] -= (p->ell_up_weight)*deltas[k];
    f_be[k] -= (p->ell_up_weight)*deltas[npts + k];
  }
  return;
}


inline void apply_up_ab_cn(const VD& deltas, const VD& old_al, const VD& old_be,
			   VD& f_al, VD& f_be, VD& cn_al, VD& cn_be, PAR *p, int npts)
{
  for (int k = 0; k < npts; ++k) {
    f_al[k] -= (p->ell_up_weight)*deltas[k];
    f_be[k] -= (p->ell_up_weight)*deltas[npts + k];
    cn_al[k] = 0.5 * (old_al[k] + f_al[k]);
    cn_be[k] = 0.5 * (old_be[k] + f_be[k]);
  }
  return;
}

inline void apply_up_abp_join(const VD& deltas, VD& f_al, VD& f_be, VD& f_ps, PAR *p, int npts)
{
  int kps = 2*npts;
  for (int k = 0; k < npts; ++k) {
    f_al[k] -= (p->ell_up_weight)*deltas[k];
    f_be[k] -= (p->ell_up_weight)*deltas[npts + k];
    f_ps[k] -= (p->ell_up_weight)*deltas[kps + k];
  }
  return;
}

inline void apply_up_abp_join_cn(const VD& deltas, const VD& old_al, const VD& old_be, const VD& old_ps,
				 VD& f_al, VD& f_be, VD& f_ps, VD& cn_al, VD& cn_be, VD& cn_ps, PAR *p, int npts)
{
  int kps = 2*npts;
  for (int k = 0; k < npts; ++k) {
    f_al[k] -= (p->ell_up_weight)*deltas[k];
    f_be[k] -= (p->ell_up_weight)*deltas[npts + k];
    f_ps[k] -= (p->ell_up_weight)*deltas[kps + k];
    cn_al[k] = 0.5 * (old_al[k] + f_al[k]);
    cn_be[k] = 0.5 * (old_be[k] + f_be[k]);
    cn_ps[k] = 0.5 * (old_ps[k] + f_ps[k]);
  }
  return;
}

#endif
