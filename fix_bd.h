/*
 *  fix_bd.h
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 1/20/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 *  -----------------------
 *  implementd brownian dynamics
 *  based on fix_nve
 *  input: needs same damping parameter as provided for the langevin fix.
 *
 */

#ifdef FIX_CLASS

FixStyle(bd,FixBD)

#else

#ifndef LMP_FIX_BD_H
#define LMP_FIX_BD_H

#include "fix.h"

namespace LAMMPS_NS {
    
    class FixBD : public Fix {
    public:
	FixBD(class LAMMPS *, int, char **);
	int setmask();
	virtual void init();
	virtual void setup(int);
	virtual void pre_force(int);
	virtual void pre_force_respa(int vflag,int,int)   { pre_force(vflag); }
	virtual void setup_pre_force(int vflag)           { pre_force(vflag); }
	virtual void setup_pre_force_respa(int vflag,int) { pre_force(vflag); }
	virtual void initial_integrate(int);
	virtual void initial_integrate_respa(int, int, int);
	virtual void final_integrate();
	virtual void final_integrate_respa(int, int);
	void reset_dt();
	
    protected:
	double dtv,dtf;
	double *step_respa;
	int mass_require;
	double t_period;  // same as in langevin
    };
    
}

#endif
#endif
