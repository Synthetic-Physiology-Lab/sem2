/*
 *  pair_semEMp.h
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 1/26/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 * ----------------------------
 * Potential for subscellular element method
 *  - based on pari_semextp
 *  - includes neighbour count based membrane detection
 * ----------------------------
 */


#ifdef PAIR_CLASS

PairStyle(semEMp,PairSemEMp)

#else

#ifndef LMP_PAIR_SEMEMP_H
#define LMP_PAIR_SEMEMP_H

#include "pair.h"

namespace LAMMPS_NS {
    
    class PairSemEMp : public Pair {
    public:
	PairSemEMp(class LAMMPS *);
	~PairSemEMp();
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
	
	double *r0,*cutDens2;
	
	double *u0,*rho,*alpha,*rhoScl,*semp1,*cut_loc;
	double **u0_inter,**rho_inter,**alpha_inter,**rhoScl_inter,**semp1_inter,**cut;
	double **offset;
	
	void allocate();
    };
    
}

#endif
#endif
