#ifndef SIMULATION_HEADER_H
#define SIMULATION_HEADER_H

#ifndef DRVAL
// doubles in MAPID r
#define NORM_FACTOR 32504
#define INDRSQ 32500
#define IN2DR 32501
#define INDR 32502
#define INDT 32503
#define RMIN 0
#define RMAX 32505
#define DRVAL 32506
#define DTVAL 32507
#define LAMVAL 32508
#define LAM2VAL 32509
#define LAM6VAL 32510
#define CSOMM 32511
#define INRMAX 32512
#define NEG2INDRSQ 32513
#define WRITEDR 32514
#define TWO_THIRDS 32515
#define TWELFTH 32516
#define FOUR_THIRDS 32517
#define FIVE_TWELFTHS 32518
#define EIGHT_PI 32519
#define TWELVE_PI 32520
#define JAC_RR 32521
#define JAC_RRM1 32522
#define JAC_RRM2 32523
#define JAC_N00 32524
#define JAC_N01 32525
#define JAC_N02 32526
#define DT_TWELVE 32527
#define CPSI_RHS 32528
#define CSOMM_RHS 32529
#define CSOMM_OLD 32530
#define DSPN_WEIGHT 32531
#define EUP_WEIGHT 32532
#define HYPTOL 32533
#define ELLTOL 32534
// integers in MAPII z
#define LASTPOINT 32535
#define SAVEPOINT 32536
#define LASTSTEP 32537
#define SAVESTEP 32538
#define NUM_POINTS 32539
#define NUM_ELL 32540
#define NUM_HYP 32541
#define RESN_FACTOR 32542
#define LASTWRITE 32543
#define WRITESHAPE 32544
#define LP_N 32545
#define LP_KL 32546
#define LP_KU 32547
#define LP_NRHS 32548
#define LP_LDAB 32549
#define LP_LDB 32550
#define MAX_ITN 32551
#endif

// libraries
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

using namespace std;

// TYPEDEFS
typedef string str;
typedef double dbl;

typedef vector<dbl> VD;
typedef map<int, double> MAPID;
typedef map<int, double> MAPII;

typedef void (*WR_FN)(const VD&, const VD&, const VD&, const VD&, const VD&,
		      const VD&, const VD&, const VD&, const VD&, const VD&,
		      const VD&, const VD&, const VD&, VD&,
		      VD&, VD&, VD&, VD&, VD&, VD&, VD&, VD&, VD&, VD&,
		      MAPID&, int, int);
typedef void (*WR_RES_FN)(const VD&, VD&, VD&, VD&, VD&, VD&, int, int);

typedef struct bbhutil_params BBHP;
typedef struct sim_fields FLDS;
typedef struct sim_writers WRS;
typedef struct sim_params PAR;


#endif
