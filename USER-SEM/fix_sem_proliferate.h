/*
 *  fix_sem_proliferate.h
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 1/6/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 *  fix to handle proliferating cells
 *  based on fix_adapt
 */

#ifdef FIX_CLASS

FixStyle(sem_proliferate,FixSemProliferate)

#else

#ifndef LMP_FIX_PROLIFERATE_H
#define LMP_FIX_PROLIFERATE_H


#include "fix.h"

namespace LAMMPS_NS {

    
class FixSemProliferate : public Fix {
public:
    
    FixSemProliferate(class LAMMPS *, int, char**);
    ~FixSemProliferate();
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
    double polyRamp, polyInc;
    double prolRamp, prolMaxInc, preprolInc, postprolInc;
    
    int maxIdCell;
    int nCell, nCreateSce, nPolymerizeSce, nSplitCell, nPrepareSplit;
    int *idCell, *ccCell, *tyCell, *nSceCell;
    bool *bCreateSce,*bSplitCell,*bPrepareSplit;
    double *dtCell, *idurCell;
    
    int polyGroupId, depolyGroupId;
    int polyGroupBitmask, depolyGroupBitmask;
    int preprolGroupId, postprolGroupId;
    int preprolGroupBitmask, postprolGroupBitmask;
    
    int sizeCell;
    
    // some private functions
    void advanceCellCycle();
    void createPolymerizeSce();
    void createSceCenter();
    void createSceRandom();
    void polymerizeSce();
    void prepareSplit();
    void splitCell();

};
}

#endif
#endif
