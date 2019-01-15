#ifndef JACOBIAN_H_INCLUDED
#define JACOBIAN_H_INCLUDED

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
#include "ekg-proc.h"


using namespace std;

//  for LAPACKE_dgbsv(): jac[ (kl + ku + 1) + (ldab - 1)*j + i ]  =  jac[ i, j ]
void set_jacCMabpslow(VD& jac, const VD& f_xi, const VD& f_pi,
		      const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, MAPID& r,
		      int npts, int kl, int ku, int ldab);
void set_jacCMabpfast(VD& jac, const VD& f_xi, const VD& f_pi,
		      const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, MAPID& r,
		      int npts, int kl, int ku, int ldab);
void set_jacCM_ab_slow(VD& jac, const VD& f_xi, const VD& f_pi,
		       const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, MAPID& r,
		       int npts, int kl, int ku, int ldab);
inline int jac_ind(int j, int k) { return (4 + j + 6*k); }
// ***********************  JACOBIAN FUNCTIONS  ***********************
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl jac_aa(const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		  PAR *p, int k)
{
  return (p->neg2indrsq) - (p->eight_pi)*sq(f_pi[k]) +
    (p->two_thirds)*pw4(f_ps[k])*sq(ddr_c(f_be,p,k) - f_be[k]*(p->r[-k])) / sq(f_al[k]);
}

inline dbl jac_aa_pm(const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, int k, int p_m)
{
  return (p->indrsq) + p_m*(p->indr)*((ddr_c(f_ps,p,k)/f_ps[k]) + (p->r[-k]));
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl jac_bb(const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, int k)
{
  return (p->neg2indrsq) - (p->r[-k])*(2*(p->r[-k]) + 6*(ddr_c(f_ps,p,k)/f_ps[k]) - (ddr_c(f_al,p,k)/f_al[k]));
}

inline dbl jac_bb_pm(const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, int k, int p_m)
{
  return (p->indrsq) + p_m*(p->in2dr)*(2*(p->r[-k]) + 6*(ddr_c(f_ps,p,k)/f_ps[k]) - (ddr_c(f_al,p,k)/f_al[k]));
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
inline dbl jac_pp(const VD& f_xi, const VD& f_pi, const VD& f_al, const VD& f_be, const VD& f_ps,
		  PAR *p, int k)
{
  return (p->neg2indrsq) + M_PI*(sq(f_xi[k]) + sq(f_pi[k])) +
    (p->five_twelfths)*pw4(f_ps[k])*sq(ddr_c(f_be,p,k) - f_be[k]*(p->r[-k])) / sq(f_al[k]);
}

inline dbl jac_pp_pm(const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p, int k, int p_m)
{
  return (p->indrsq) + p_m*(p->indr)*(p->r[-k]);
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// **********************************************************
// **********************************************************
//                    POPULATING JACOBIAN
// **********************************************************
// **********************************************************
////////////////////////////////////////////////////////////////////////////////////////////////

void set_jacCMabpfast(VD& jac, const VD& f_xi, const VD& f_pi,
		      const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p,
		      MAPID& r, int npts, int kl, int ku, int ldab)
{
  int one_past_last = npts - 1;
  int j = 0, jbe = npts, jps = 2*npts;
  // ROW 0, COL 0
  jac[jac_ind(j,j)] = (p->jacN00);
  jac[jac_ind(jbe,jbe)] = 1;
  jac[jac_ind(jps,jps)] = (p->jacN00);
  // ROW 0, COL 1
  jac[jac_ind(j,j + 1)] = (p->jacN01);
  jac[jac_ind(jbe,jbe + 1)] = 0;
  jac[jac_ind(jps,jps + 1)] = (p->jacN01);
  // ROW 0, COL 2
  jac[jac_ind(j,j + 2)] = (p->jacN02);
  jac[jac_ind(jbe,jbe + 2)] = 0;
  jac[jac_ind(jps,jps + 2)] = (p->jacN02);

  double drin_p, drin_m, dps_ps, dlogp6_a, p4db2_a2;
  for (j = 1; j < one_past_last; ++j) {
    jbe = j + npts;
    jps = jbe + npts;

    drin_p = (p->indrsq) + r[-j]*(p->indr);
    drin_m = (p->indrsq) - r[-j]*(p->indr);
    dps_ps = ddr_c(f_ps,p,j) / f_ps[j];
    dlogp6_a = 6*dps_ps - (ddr_c(f_al,p,j)/f_al[j]);
    p4db2_a2 = sq(f_ps[j]) * (ddr_c(f_be,p,j) - r[-j]*f_be[j]) / f_al[j];
    p4db2_a2 = sq(p4db2_a2);    
    
    // ROW j, COL j-1
    jac[jac_ind(j,j - 1)] = drin_m - r[-j]*dps_ps;
    jac[jac_ind(jbe,jbe - 1)] = drin_m - (p->in2dr)*dlogp6_a;
    jac[jac_ind(jps,jps - 1)] = drin_m;
    // ROW j, COL j
    jac[jac_ind(j,j)] = (p->neg2indrsq) + (p->two_thirds)*p4db2_a2 - (p->eight_pi)*sq(f_pi[j]);
    jac[jac_ind(jbe,jbe)] = (p->neg2indrsq) - r[-j]*(2*r[-j] + dlogp6_a);
    jac[jac_ind(jps,jps)] = (p->neg2indrsq) + (p->five_twelfths)*p4db2_a2 + M_PI*(sq(f_xi[j]) + sq(f_pi[j]));
    // ROW j+1, COL j
    jac[jac_ind(j,j + 1)] = drin_p + r[-j]*dps_ps;
    jac[jac_ind(jbe,jbe + 1)] = drin_p + (p->in2dr)*dlogp6_a;
    jac[jac_ind(jps,jps + 1)] = drin_p;
  }
  
  j = one_past_last;
  jbe = j + npts;
  jps = jbe + npts;
  // ROW lastpt, COL lastpt-2
  jac[jac_ind(j,j - 2)] = (p->jacRRm2);
  jac[jac_ind(jbe,jbe - 2)] = (p->jacRRm2);
  jac[jac_ind(jps,jps - 2)] = (p->jacRRm2);
  // ROW lastpt, COL lastpt-1
  jac[jac_ind(j,j - 1)] = (p->jacRRm1);
  jac[jac_ind(jbe,jbe - 1)] = (p->jacRRm1);
  jac[jac_ind(jps,jps - 1)] = (p->jacRRm1);
  // ROW lastpt, COL lastpt
  jac[jac_ind(j,j)] = (p->jacRR);
  jac[jac_ind(jbe,jbe)] = (p->jacRR);
  jac[jac_ind(jps,jps)] = (p->jacRR);
  return;
}

void set_jacCMabpslow(VD& jac, const VD& f_xi, const VD& f_pi,
		      const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p,
		      MAPID& r, int npts, int kl, int ku, int ldab)
{
  int one_past_last = npts - 1;
  int j = 0, jbe = npts, jps = 2*npts;
  // ROW 0, COL 0
  jac[jac_ind(j,j)] = (p->jacN00);
  jac[jac_ind(jbe,jbe)] = 1;
  jac[jac_ind(jps,jps)] = (p->jacN00);
  // ROW 0, COL 1
  jac[jac_ind(j,j + 1)] = (p->jacN01);
  jac[jac_ind(jbe,jbe + 1)] = 0;
  jac[jac_ind(jps,jps + 1)] = (p->jacN01);
  // ROW 0, COL 2
  jac[jac_ind(j,j + 2)] = (p->jacN02);
  jac[jac_ind(jbe,jbe + 2)] = 0;
  jac[jac_ind(jps,jps + 2)] = (p->jacN02);
  
  for (j = 1; j < one_past_last; ++j) {
    jbe = j + npts;
    jps = jbe + npts;
    // ROW j, COL j-1
    jac[jac_ind(j,j - 1)] = jac_aa_pm(f_al, f_be, f_ps, p, j, -1);
    jac[jac_ind(jbe,jbe - 1)] = jac_bb_pm(f_al, f_be, f_ps, p, j, -1);
    jac[jac_ind(jps,jps - 1)] = jac_pp_pm(f_al, f_be, f_ps, p, j, -1);
    // ROW j, COL j
    jac[jac_ind(j,j)] = jac_aa(f_xi, f_pi, f_al, f_be, f_ps, p, j);
    jac[jac_ind(jbe,jbe)] = jac_bb(f_al, f_be, f_ps, p, j);
    jac[jac_ind(jps,jps)] = jac_pp(f_xi, f_pi, f_al, f_be, f_ps, p, j);
    // ROW j+1, COL j
    jac[jac_ind(j,j + 1)] = jac_aa_pm(f_al, f_be, f_ps, p, j, 1);
    jac[jac_ind(jbe,jbe + 1)] = jac_bb_pm(f_al, f_be, f_ps, p, j, 1);
    jac[jac_ind(jps,jps + 1)] = jac_pp_pm(f_al, f_be, f_ps, p, j, 1);
  }
  j = one_past_last;
  jbe = j + npts;
  jps = jbe + npts;
  // ROW lastpt, COL lastpt-2
  jac[jac_ind(j,j - 2)] = (p->jacRRm2);
  jac[jac_ind(jbe,jbe - 2)] = (p->jacRRm2);
  jac[jac_ind(jps,jps - 2)] = (p->jacRRm2);
  // ROW lastpt, COL lastpt-1
  jac[jac_ind(j,j - 1)] = (p->jacRRm1);
  jac[jac_ind(jbe,jbe - 1)] = (p->jacRRm1);
  jac[jac_ind(jps,jps - 1)] = (p->jacRRm1);
  // ROW lastpt, COL lastpt
  jac[jac_ind(j,j)] = (p->jacRR);
  jac[jac_ind(jbe,jbe)] = (p->jacRR);
  jac[jac_ind(jps,jps)] = (p->jacRR);
  return;
}

void set_jacCM_ab_slow(VD& jac, const VD& f_xi, const VD& f_pi,
		       const VD& f_al, const VD& f_be, const VD& f_ps, PAR *p,
		       MAPID& r, int npts, int kl, int ku, int ldab)
{
  int one_past_last = npts - 1;
  int j = 0, jbe = npts;
  // ROW 0, COL 0
  jac[jac_ind(j,j)] = (p->jacN00);
  jac[jac_ind(jbe,jbe)] = 1;
  // ROW 0, COL 1
  jac[jac_ind(j,j + 1)] = (p->jacN01);
  jac[jac_ind(jbe,jbe + 1)] = 0;
  // ROW 0, COL 2
  jac[jac_ind(j,j + 2)] = (p->jacN02);
  jac[jac_ind(jbe,jbe + 2)] = 0;
  for (j = 1; j < one_past_last; ++j) {
    jbe = j + npts;
    // ROW j, COL j-1
    jac[jac_ind(j,j - 1)] = jac_aa_pm(f_al, f_be, f_ps, p, j, -1);
    jac[jac_ind(jbe,jbe - 1)] = jac_bb_pm(f_al, f_be, f_ps, p, j, -1);
    // ROW j, COL j
    jac[jac_ind(j,j)] = jac_aa(f_xi, f_pi, f_al, f_be, f_ps, p, j);
    jac[jac_ind(jbe,jbe)] = jac_bb(f_al, f_be, f_ps, p, j);
    // ROW j+1, COL j
    jac[jac_ind(j,j + 1)] = jac_aa_pm(f_al, f_be, f_ps, p, j, 1);
    jac[jac_ind(jbe,jbe + 1)] = jac_bb_pm(f_al, f_be, f_ps, p, j, 1);
  }
  j = one_past_last;
  jbe = j + npts;
  // ROW lastpt, COL lastpt-2
  jac[jac_ind(j,j - 2)] = (p->jacRRm2);
  jac[jac_ind(jbe,jbe - 2)] = (p->jacRRm2);
  // ROW lastpt, COL lastpt-1
  jac[jac_ind(j,j - 1)] = (p->jacRRm1);
  jac[jac_ind(jbe,jbe - 1)] = (p->jacRRm1);
  // ROW lastpt, COL lastpt
  jac[jac_ind(j,j)] = (p->jacRR);
  jac[jac_ind(jbe,jbe)] = (p->jacRR);
  return;
}




#endif
