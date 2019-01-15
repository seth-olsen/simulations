#ifndef FDA_FNS_H_INCLUDED
#define FDA_FNS_H_INCLUDED

#include <iostream>
#include <vector> // for everything
#include <cmath> // for ICs
#include <map>
#include <string>
#include "sim-structs.h"
#include "ekg-header.h"

using namespace std;

// x^n and 1/x^n
inline dbl sq(dbl x) { return (x*x); }
inline dbl pw3(dbl x) { return (x*x*x); }
inline dbl pw4(dbl x) { return (x*x*x*x); }
inline dbl pw5(dbl x) { return (x*x*x*x*x); }
inline dbl pw6(dbl x) { return (x*x*x*x*x*x); }

inline dbl norm_inf(const VD& vec)
{
  return max( *max_element(vec.begin(), vec.end()),
	     -(*min_element(vec.begin(), vec.end())) );
}

// ******** DIFFERENCES ********

inline dbl d_c(const VD& u, int ind)
{ return u[ind+1] - u[ind-1]; }

inline dbl d_f(const VD& u, int ind)
{ return -3*u[ind] + 4*u[ind+1] - u[ind+2]; }

inline dbl d_b(const VD& u, int ind)
{ return 3*u[ind] - 4*u[ind-1] + u[ind-2]; }

inline dbl d2_c(const VD& u, int ind)
{ return u[ind-1] - 2*u[ind] + u[ind+1]; }

inline dbl d2_f(const VD& u, int ind)
{ return 2*u[ind] - 5*u[ind+1] + 4*u[ind+2] - u[ind+3]; }

inline dbl d2_b(const VD& u, int ind)
{ return 2*u[ind] - 5*u[ind-1] + 4*u[ind-2] - u[ind-3]; }

inline dbl dlog_c(const VD& u, int ind)
{ return log(u[ind+1] / u[ind-1]); }
inline dbl dlog_f(const VD& u, int ind)
{ return -3*log(u[ind]) + 4*log(u[ind+1]) - log(u[ind+2]); }
inline dbl dlog_b(const VD& u, int ind)
{ return 3*log(u[ind]) - 4*log(u[ind+1]) + log(u[ind+2]); }

// DERIVATIVES

inline dbl ddr_c(const VD& u, PAR *p, int ind)
{ return (p->in2dr)*(u[ind+1] - u[ind-1]); }

inline dbl ddr_f(const VD& u, PAR *p, int ind)
{ return (p->in2dr)*(-3*u[ind] + 4*u[ind+1] - u[ind+2]); }

inline dbl ddr_b(const VD& u, PAR *p, int ind)
{ return (p->in2dr)*(3*u[ind] - 4*u[ind-1] + u[ind-2]); }

inline dbl ddr2_c(const VD& u, PAR *p, int ind)
{ return (p->indrsq)*(u[ind-1] - 2*u[ind] + u[ind+1]); }

inline dbl ddr2_f(const VD& u, PAR *p, int ind)
{ return (p->indrsq)*(2*u[ind] - 5*u[ind+1] + 4*u[ind+2] - u[ind+3]); }

inline dbl ddr2_b(const VD& u, PAR *p, int ind)
{ return (p->indrsq)*(2*u[ind] - 5*u[ind-1] + 4*u[ind-2] - u[ind-3]); }

inline dbl ddrlog_c(const VD& u, PAR *p, int ind)
{ return (p->in2dr)*log(u[ind+1] / u[ind-1]); }
inline dbl ddrlog_f(const VD& u, PAR *p, int ind)
{ return (p->in2dr)*(-3*log(u[ind]) + 4*log(u[ind+1]) - log(u[ind+2])); }
inline dbl ddrlog_b(const VD& u, PAR *p, int ind)
{ return (p->in2dr)*(3*log(u[ind]) - 4*log(u[ind+1]) + log(u[ind+2])); }

// d(r*u)/dr * 2*dr
inline dbl d_ru_c(const VD& u, int ind, dbl dr, dbl r)
{ return u[ind+1]*(r+dr) - u[ind-1]*(r-dr); }

inline dbl d_ru_f(const VD& u, int ind, dbl dr, dbl r)
{ return -3*u[ind]*r + 4*u[ind+1]*(r+dr) - u[ind+2]*(r+2*dr); }

inline dbl d_ru_b(const VD& u, int ind, dbl dr, dbl r)
{ return 3*u[ind]*r - 4*u[ind-1]*(r-dr) + u[ind-2]*(r-2*dr); }

// d(u/r)/dr * 2*dr
inline dbl d_urinv_c(const VD& u, int ind, dbl dr, dbl r)
{ return u[ind+1]/(r+dr) - u[ind-1]/(r-dr); }

inline dbl d_urinv_f(const VD& u, int ind, dbl dr, dbl r)
{ return -3*u[ind]/r + 4*u[ind+1]/(r+dr) - u[ind+2]/(r+2*dr); }

inline dbl d_urinv_b(const VD& u, int ind, dbl dr, dbl r)
{ return 3*u[ind]/r - 4*u[ind-1]/(r-dr) + u[ind-2]/(r-2*dr); }

// d(r^2*u)/dr * 2*dr
inline dbl d_r2u_c(const VD& u, int ind, dbl dr, dbl r)
{ return u[ind+1]*sq(r+dr) - u[ind-1]*sq(r-dr); }

inline dbl d_r2u_f(const VD& u, int ind, dbl dr, dbl r)
{ return -3*u[ind]*sq(r) + 4*u[ind+1]*sq(r+dr) - u[ind+2]*sq(r+2*dr); }

inline dbl d_r2u_b(const VD& u, int ind, dbl dr, dbl r)
{ return 3*u[ind]*sq(r) - 4*u[ind-1]*sq(r-dr) + u[ind-2]*sq(r-2*dr); }

// d(u/r^2)/dr * 2*dr
inline dbl d_ur2inv_c(const VD& u, int ind, dbl dr, dbl r)
{ return u[ind+1]/sq(r+dr) - u[ind-1]/sq(r-dr); }

inline dbl d_ur2inv_f(const VD& u, int ind, dbl dr, dbl r)
{ return -3*u[ind]/sq(r) + 4*u[ind+1]/sq(r+dr) - u[ind+2]/sq(r+2*dr); }

inline dbl d_ur2inv_b(const VD& u, int ind, dbl dr, dbl r)
{ return 3*u[ind]/sq(r) - 4*u[ind-1]/sq(r-dr) + u[ind-2]/sq(r-2*dr); }

// d(r^2*u^4)/dr * 2*dr
inline dbl d_r2u4_c(const VD& u, int ind, dbl dr, dbl r)
{ return ( sq(r+dr)*pw4(u[ind+1]) - sq(r-dr)*pw4(u[ind-1]) ); }

inline dbl d_r2u4_f(const VD& u, int ind, dbl dr, dbl r)
{ return ( -3*sq(r)*pw4(u[ind]) + 4*sq(r+dr)*pw4(u[ind+1]) - sq(r+2*dr)*pw4(u[ind+2]) ); }

inline dbl d_r2u4_b(const VD& u, int ind, dbl dr, dbl r)
{ return ( 3*sq(r)*pw4(u[ind]) - 4*sq(r-dr)*pw4(u[ind-1]) + sq(r-2*dr)*pw4(u[ind-2]) ); }

// d(u/v^2)/dr * 2*dr
inline dbl d_uv2inv_c(const VD& u, const VD& v, int ind)
{ return ( u[ind+1]/sq(v[ind+1]) - u[ind-1]/sq(v[ind-1]) ); }

inline dbl d_uv2inv_f(const VD& u, const VD& v, int ind)
{ return ( -3*u[ind]/sq(v[ind]) + 4*u[ind+1]/sq(v[ind+1]) - u[ind+2]/sq(v[ind+2]) ); }

inline dbl d_uv2inv_b(const VD& u, const VD& v, int ind)
{ return ( 3*u[ind]/sq(v[ind]) - 4*u[ind-1]/sq(v[ind-1]) + u[ind-2]/sq(v[ind-2]) ); }


// CRANK-NICHOLSON
inline dbl cn(const VD& old_f, const VD& f, int ind)
{
  return 0.5*(old_f[ind] + f[ind]);
}

inline void set2_cn(const VD& old_f1, const VD& old_f2,
		    const VD& f1, const VD& f2,
		    VD& cn_f1, VD& cn_f2, int npts)
{
  for (int k = 0; k < npts; ++k) {
      cn_f1[k] = 0.5*(old_f1[k] + f1[k]);
      cn_f2[k] = 0.5*(old_f2[k] + f2[k]);
  }
  return;
}
inline void set3_cn(const VD& old_f1, const VD& old_f2, const VD& old_f3,
		    const VD& f1, const VD& f2, const VD& f3,
		    VD& cn_f1, VD& cn_f2, VD& cn_f3, int npts)
{
  for (int k = 0; k < npts; ++k) {
      cn_f1[k] = 0.5*(old_f1[k] + f1[k]);
      cn_f2[k] = 0.5*(old_f2[k] + f2[k]);
      cn_f3[k] = 0.5*(old_f3[k] + f3[k]);
  }
  return;
}

// *********** BOUNDARY CONDITIONS ****************


inline void dirichlet0(VD& field)
{ field[0] = 0; }

inline void neumann0(VD& field)
{ field[0] = (4*field[1] - field[2]) / 3.0; }

inline void sommerfeld(const VD& oldfield, VD& field, PAR *p,
		       int ind)
{
  field[ind] = (p->csomm_rhs)*( (p->lam)*(field[ind-1] + oldfield[ind-1])
			      - 0.25*(p->lam)*(field[ind-2] + oldfield[ind-2])
			      + (p->csomm_old)*oldfield[ind] );
  return;
}

inline dbl dirichlet0res(const VD& field)
{ return field[0]; }

inline dbl neumann0res(const VD& field, PAR *p)
{ return ddr_f(field, p, 0); }

inline dbl sommerfeldres(const VD& oldfield, const VD& field,
			 PAR *p, int ind)
{
  return 0.5*(field[ind] + oldfield[ind]) +
    (p->r[ind])*( (p->indt)*(field[ind] - oldfield[ind])
	    + 0.5*(ddr_b(field,p,ind) + ddr_b(oldfield,p,ind)) );
}


// ********** DISSIPATION FUNCTIONS **************


// kreiss-oliger dissipation (p.23 choptuik notes)
inline dbl dissipate(dbl eps, const VD& u, int ind)
{ return -0.0625 * eps * ( u[ind-2] - 4*u[ind-1] + 6*u[ind]
			       - 4*u[ind+1] + u[ind+2] ); }

inline dbl symdiss1(dbl eps, const VD& u)
{ return -0.0625 * eps * ( u[3] - 4*u[2] + 7*u[1] - 4*u[0] ); }

inline dbl antidiss1(dbl eps, const VD& u)
{ return -0.0625 * eps * ( u[3] - 4*u[2] + 5*u[1] - 4*u[0] ); }

inline dbl symdiss0(dbl eps, const VD& u)
{ return -0.0625 * eps * ( 2*u[2] - 8*u[1] + 6*u[0] ); }

inline dbl antidiss0(dbl eps, const VD& u)
{ return -0.0625 * eps * ( 6*u[0] ); }



#endif

