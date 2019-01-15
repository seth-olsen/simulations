#ifndef EKG_PROC_H_INCLUDED
#define EKG_PROC_H_INCLUDED

#ifndef INDRSQ
#define INDRSQ 32500
#define IN2DR 32501
#define INDT 32502
#define RMIN 32503
#define RMAX 32504
#define DRVAL 32505
#define DTVAL 32506
#define LAMVAL 32507
#define LAM2VAL 32508
#define LAM6VAL 32509
#define CSOMM 32510
#define INRMAX 32511
#endif

#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <vector> // for everything
#include <cmath> // for ICs
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf
#include "fda-io.h"
#include "fda-fns.h"
#include "ekg-fns.h"
#include "jacobian.h"
#include "ekg-clean.h"
#include "lapacke.h"

using namespace std;

// dissipation functions
void dissipationNB_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi,
		      int one_past_last, dbl dspn);
void dissipationB_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi,
		     int one_past_last, dbl dspn);
void dissipationNB2_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi,
		       int one_past_last, dbl dspn);
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// at ind next to boundaries can call dissipate on ind+/-1 or ind+/-2
//************ CHECKED: Mon. 10/22 ************
void dissipationNB_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi,
		      int one_past_last, dbl dspn)
{
  f_xi[1] += antidiss1(dspn, old_xi);
  f_pi[1] += symdiss1(dspn, old_pi);
  for (int k = 2; k < one_past_last; ++k) {
    f_xi[k] += dissipate(dspn, old_xi, k);
    f_pi[k] += dissipate(dspn, old_pi, k);
  }
  return;
}
void dissipationNB2_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi,
		       int one_past_last, dbl dspn)
{
  for (int k = 2; k < one_past_last; ++k) {
    f_xi[k] += dissipate(dspn, old_xi, k);
    f_pi[k] += dissipate(dspn, old_pi, k);
  }
  return;
}
//************ CHECKED: Mon. 10/22 ************
void dissipationB_xp(const VD& old_xi, const VD& old_pi, VD& f_xi, VD& f_pi,
		     int one_past_last, dbl dspn)
{
  f_pi[0] += symdiss0(dspn, old_pi);
  f_xi[1] += antidiss1(dspn, old_xi);
  f_pi[1] += symdiss1(dspn, old_pi);
  for (int k = 2; k < one_past_last; ++k) {
    f_xi[k] += dissipate(dspn, old_xi, k);
    f_pi[k] += dissipate(dspn, old_pi, k);
  }
  return;
}
  
#endif

