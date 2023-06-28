
/*
 *  fix_sem_PMN.h
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 1/26/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 * ++++++++++++++++++++++++++++++++++++++++++++++
 * - fix to handle proliferation and migration with a nucleus
 * - based on fix_sem_PM
 * ++++++++++++++++++++++++++++++++++++++++++++++
 */



#ifdef FIX_CLASS

FixStyle(sem_PMN,FixSemPMN)

#else

#ifndef LMP_FIX_SEMPMN_H
#define LMP_FIX_SEMPMN_H


#include "fix.h"

namespace LAMMPS_NS {
    
    
    class FixSemPMN : public Fix {
    public:
	
	FixSemPMN(class LAMMPS *, int, char**);
	~FixSemPMN();
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
	double polyNucRamp, polyNucInc;
	
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
	// for splitting nucleus
	int polyNucGroupId, polyNucGroupBitmask;
	
	int sizeCell;
	
	// MIGRATION
	double ** polarity;
	
	// Cell Cycle
	void advanceCell();
	
	// Nucleus centering
	void nuc_centre();
	
	// Growth
	void createSce();
	void createSceCenter();
	void createSceRandom();
	void prepareSplit();
	void prepareSplitInflateElements();
	void prepareSplitNuc();
	void splitCell();
	void splitCellGeometric();
	void splitCellNuc();
	
	//Polymerization
	void polymerizeSce();
	
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
