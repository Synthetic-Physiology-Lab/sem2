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

PairStyle(semp,PairSemp)

#else

#ifndef LMP_PAIR_SEMP_H
#define LMP_PAIR_SEMP_H

#include "pair.h"

namespace LAMMPS_NS {
	
	class PairSemp : public Pair {
	public:
		PairSemp(class LAMMPS *);
		~PairSemp();
		void compute(int, int);
		void settings(int, char **);
		void coeff(int, char **);
		double init_one(int, int);
		void write_restart(FILE *);
		void read_restart(FILE *);
		void write_restart_settings(FILE *);
		void read_restart_settings(FILE *);
		double single(int, int, int, int, int, int, double, double, double, double &);
		
	protected:
		double cut_global;
		double **cut;
		
		double **r0;
		double *cut_loc,*r02,*rhoScl,*u0,*semp1;
		double *rho_inter;
		double **u0_inter,**r02_inter,**rhoScl_inter,**semp1_inter;
		double **offset;
		
		void allocate();
	};
	
}

#endif
#endif
