/*
 *  fix_sem_polymerize.h
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 3/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 *  fix to handle prolymerization of atom types
 *  based on fix_polymerize
 */

#ifdef FIX_CLASS

FixStyle(sem_polymerize,FixSemPolymerize)

#else

#ifndef LMP_FIX_POLYMERIZE_H
#define LMP_FIX_POLYMERIZE_H


#include "fix.h"

namespace LAMMPS_NS {
    
    
    class FixSemPolymerize : public Fix {
    public:
	
	FixSemPolymerize(class LAMMPS *, int, char**);
	~FixSemPolymerize();
	int setmask();		// determines when fix is called during the time step
	void init();		// initialization before run
	
	void final_integrate();
	void end_of_step();
	
	void pre_force(int);        // called before forces are being calculated
	// -ramp up polymerization constant (maybe not here)
	// -split cell
	
	void post_run();		// clean up after run
	
    private:
	// some private variables
	double polyDur, polyInc;
	
	int groupid;
	
	int nPolymerizeSce;
	
	int polyGroupId;
	int polyGroupBitmask;
	
	// some private functions
	void polymerizeSce();	
    };
}

#endif
#endif
