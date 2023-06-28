/*
 *  fix_sem_polymerize.cpp
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 3/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_sem_polymerize.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "kspace.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

#include "domain.h"
#include "group.h"
#include "atom_vec_sem.h"
#include "eig3.h"

// define __SEM_DEBUG__ to get additional printfs
//#define __SEM_DEBUG__

using namespace LAMMPS_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
 input: fix $fixid all sem_proliferate $Nevery $dur
 for now.. maybe list of cell ids
 ---------------------------------------------------------------------- */
FixSemPolymerize::FixSemPolymerize(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
    
    if (narg < 5) error->all(FLERR,"Illegal fix sem_proliferate command - not enough parameters");
    nevery = utils::inumeric(FLERR,arg[3],false,lmp);
    if (nevery < 0) error->all(FLERR,"Illegal fix sem_proliferate command - nevery negative");
    
    dynamic_group_allow = 1;
    
    // cell cycle times
    groupid = group->find(arg[1]);
    polyDur = atof(arg[4]);
        
    // enable force reneighboring
    force_reneighbor = 1;
    next_reneighbor = -1;
}
/* ----------------------------------------------------------------------
 Destruct
 ---------------------------------------------------------------------- */
FixSemPolymerize::~FixSemPolymerize()
{
}
/* ----------------------------------------------------------------------
 When to invoce fix during execution pipeline
 ---------------------------------------------------------------------- */
int FixSemPolymerize::setmask()
{
    int mask = 0;
    mask |= END_OF_STEP;
    return mask;
}
/* ----------------------------------------------------------------------
 Init things:
 ---------------------------------------------------------------------- */
void FixSemPolymerize::init()
{
    // check some things
    if (!atom->sem_flag) error->all(FLERR,"Fix adapt requires atom of sem style");
    
    polyInc=update->dt*((double)nevery)/polyDur;
	     
    int nPolymerizeSceLoc = 0;
    double *p = atom->p;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int *polyflag = (int *) memory->smalloc(nlocal*sizeof(int),"fix/semPM:flag");

    // init atoms
    int groupMask=group->bitmask[groupid];
    
    for(int i=0;i<nlocal;i++){
	polyflag[i]=0;
	if(mask[i]==groupMask && p[i] < 1){
	    polyflag[i]=1;
	    nPolymerizeSceLoc++;
	}
    }
    
    // distribute info
    MPI_Allreduce(&nPolymerizeSceLoc,&nPolymerizeSce,1,MPI_INT,MPI_SUM,world);
    
    //create groups for poly/depoly
    group->create("poly",polyflag);
    polyGroupId=group->find("poly");
    polyGroupBitmask=group->bitmask[polyGroupId];
        
    memory->sfree(polyflag);
    
}
/* ---------------------------------------------------------------------
 - bla
 ---------------------------------------------------------------------- */
void FixSemPolymerize::final_integrate(){
    
}

/* ----------------------------------------------------------------------
 What to do after step
 ---------------------------------------------------------------------- */
void FixSemPolymerize::end_of_step(){
    if (nevery == 0) return;
    if (update->ntimestep % nevery) return;
    if(nPolymerizeSce>0){
	polymerizeSce();
    }
}


/* ----------------------------------------------------------------------
 What to do before forces are calculated
 ---------------------------------------------------------------------- */
void FixSemPolymerize::pre_force(int){
    
}
/* ----------------------------------------------------------------------
 What to do after run
 ---------------------------------------------------------------------- */
void FixSemPolymerize::post_run(){
    
}
/* ----------------------------------------------------------------------
 polymerize in group
 ---------------------------------------------------------------------- */
void FixSemPolymerize::polymerizeSce(){
#ifdef __SEM_DEBUG__
    printf(" - Polymerizing %d SCEs \n",nPolymerizeSce);
#endif
    int nPolymerizeSceLoc=0;
    int nlocal = atom->nlocal;
    int *mask  = atom->mask;
    double *p = atom->p;
    
    int me;
    MPI_Comm_rank(world,&me);
    
    for (int i=0;i<nlocal;i++){
	//	printf("--%d | proc: %d, particle %d, location: %f %f %f\n",update->ntimestep, me, atom->tag[i],x[i][0],x[i][1],x[i][2]);
	// use loop to update polymerization factors
	if(mask[i] & polyGroupBitmask){
	    nPolymerizeSceLoc++;
	    p[i]+=polyInc;
	    if(p[i]>1.0){
		p[i]=1.0;
		mask[i]^=polyGroupBitmask;
		nPolymerizeSceLoc--;
	    }
	}
    }
    // update atom->natoms;
    MPI_Allreduce(&nPolymerizeSceLoc,&nPolymerizeSce,1,MPI_INT,MPI_SUM,world);
}
