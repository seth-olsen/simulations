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
#include "sim-init.h"

int main(int argc, char **argv)
{
  // INITITALIZATION
  PAR p;
  FLDS f;
  WRS wr;
  int err_code = params_init(&p, argc, argv);
  if (err_code != 0) {
    cout << "\nPARAM INIT error code = " << error_code << endl;
    return error_code;
  }
  err_code = fields_init(&f, &p);
  if (err_code != 0) {
    cout << "\nFIELD INIT error code = " << error_code << endl;
    return error_code;
  }
  vector<BBHP *> writer_vec = writers_init(&wr, &f, &p);
  gft_set_multi();
  // DO SOME CHECKS HERE
  for (int i = 0; i < p.nsteps; ++i) {
    if (i % p.save_step == 0) {
      write_bbhp_vec(writer_vec, &p);
      if (i % p->check_diagnostics == 0) {
	write_diagnostics(&wr, &f, &p);
      }
    }
    // SOLVE FOR NEXT STEP
    fields_step(&f, &p);
  }
  // WRITE LAST STEP
  if (i % p.save_step == 0) {
    write_bbhp_vec(writer_vec, &p);
    if (i % p->check_diagnostics == 0) {
      write_diagnostics(&wr, &f, &p);
    }
  }
  gft_close_all();
  return 0;
}
    
  
