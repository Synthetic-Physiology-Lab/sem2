/*
 *  pari_semp.h
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 12/9/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 * ----------------------------
 * Potential for subscellular element method
 * ----------------------------
 */

#ifdef PAIR_CLASS

PairStyle(semextp,PairSemextp)

#else

#ifndef LMP_PAIR_SEMEXTP_H
#define LMP_PAIR_SEMEXTP_H

#include "pair.h"

namespace LAMMPS_NS {
 
class PairSemextp : public Pair {
public:
  PairSemextp(class LAMMPS *);
  ~PairSemextp();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, int, int, double, double, double, double, double &);
  
protected:
  double cut_global;
  
  double *r0;
  
  double *u0,*rho,*alpha,*rhoScl,*semp1,*cut_loc;
  double **u0_inter,**rho_inter,**alpha_inter,**rhoScl_inter,**semp1_inter,**cut;
  double **offset;
  
  void allocate();
};
 
}

#endif
#endif
