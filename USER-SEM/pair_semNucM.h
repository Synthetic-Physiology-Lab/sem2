/*
 *  pair_semNucM.h
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 1/26/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 * ----------------------------
 * Potential for subscellular element method
 *  - MUST: nucleus elements as atom types 2
 *  - Morse potential
 *  - potential acting inbetween nuclei and SCE's of same cell
 *  - based on pari_semEMp
 *  - includes neighbour count based membrane detection??
 * ----------------------------
 */


#ifdef PAIR_CLASS

PairStyle(semNucM,PairSemNucM)

#else

#ifndef LMP_PAIR_SEMNUCM_H
#define LMP_PAIR_SEMNUCM_H

#include "pair.h"

namespace LAMMPS_NS {
    
    class PairSemNucM : public Pair {
    public:
	PairSemNucM(class LAMMPS *);
	~PairSemNucM();
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
	
	double **r0,**cutDens2;
	
	double **u0,**rho,**alpha,**rShift,**cut;
	double **rhoScl,**semp1,**rShiftSq,**rMinSq;  // derrived parameters
	double **offset;
	
	void allocate();
    };
    
}

#endif
#endif
