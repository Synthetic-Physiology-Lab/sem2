/*
 *  fix_sem_PM.h
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 1/26/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 * ++++++++++++++++++++++++++++++++++++++++++++++
 * - fix to handle proliferation and migration
 * - based on fix_proliferation
 * ++++++++++++++++++++++++++++++++++++++++++++++
 */



#ifdef FIX_CLASS

FixStyle(sem_PM,FixSemPM)

#else

#ifndef LMP_FIX_SEMPM_H
#define LMP_FIX_SEMPM_H


#include "fix.h"

namespace LAMMPS_NS {
    
    
    class FixSemPM : public Fix {
    public:
	
	FixSemPM(class LAMMPS *, int, char**);
	~FixSemPM();
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
	double CC_dur[4];
	double CP_dur[3];
	double polyPRamp, polyPInc;
	double polyMRamp, polyMInc;
	double prolRamp, prolMaxInc, preprolInc, postprolInc;
	
	int maxIdCell;
	int nCell, nCreateSce, nPolymerizeSce, nSplitCell, nPrepareSplit, nPrepareMig, nMig;

	int *idCell, *ccCell, *cpCell, *tyCell, *nSceCell;
	bool *bCreateSce,*bSplitCell,*bPrepareSplit,*bPrepareMig,*bMig;
	double *dtCell, *idurCell;
	
	// masks
	// for proliferation
	int polyPGroupId, depolyPGroupId;
	int polyPGroupBitmask, depolyPGroupBitmask;
	// for migration
	int polyMGroupId, depolyMGroupId;
	int polyMGroupBitmask, depolyMGroupBitmask;
	// for splitting
	int preprolGroupId, postprolGroupId;
	int preprolGroupBitmask, postprolGroupBitmask;
	
	int sizeCell;
	
	// MIGRATION
	double ** polarity;
	
	// some private functions
	void advanceCell();
	void createSce();
	void createSceCenter();
	void createSceRandom();
	void polymerizeSce();
	void prepareSplit();
	void splitCell();
	
	void tagMembrane();
	void tagMembraneDens();
	void updatePolarity();
	
	void prepareMig();
	void migrate();
	void migratePolarity();
	void deleteSce();
	
	// HELPERS
	void updatePolarity(int,int *,int *,double **,double *);
	void cellCenter(int,bool *,double **, int *, int *);
    };
}

#endif
#endif
