
/*
 *  fix_sem_PM.cpp
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 1/6/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 * ++++++++++++++++++++++++++++++++++++++++++++++
 * - fix to handle proliferation and migration
 * - based on fix_proliferation
 * =================================
 * - grow registered cells to proper size
 * - equilibrate for a bit
 * - switch to migration
 *
 *  Edited by Sandipan on 28/06/2023
 * ++++++++++++++++++++++++++++++++++++++++++++++
 */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_sem_PM.h"
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

#define __SEM_DEBUG__NOT

#define __PROL__NOT

using namespace LAMMPS_NS;
using namespace FixConst;

enum CP {CP_Q, CP_M, CP_P};
enum CC {CC_G0,CC_I,CC_M1,CC_M2};


struct mat33{
    double v[3][3];
};
/* ----------------------------------------------------------------------
 input: fix $fixid all sem_proliferate $Nevery $nSceperCell $tG0 $tGI $tM1 $cellID+
 for now.. maybe list of cell ids
 ---------------------------------------------------------------------- */
FixSemPM::FixSemPM(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
    int fixargs=14;
    printf("narg: %d\n",narg);
    if (narg < fixargs+1) error->all(FLERR,"Illegal fix sem_PM command - not enough parameters");
    nevery = utils::inumeric(FLERR,arg[3],false,lmp);
    if (nevery < 0) error->all(FLERR,"Illegal fix sem_PM command - nevery negative");
    
    dynamic_group_allow = 1;
    
    // cell cycle times
    sizeCell = atoi(arg[4]);
    CP_dur[CP_P] = atof(arg[5]);
    CC_dur[CC_G0]=atof(arg[6]);
    CC_dur[CC_I]=atof(arg[7]);
    CC_dur[CC_M1]=atof(arg[8]);
    CC_dur[CC_M2]=atof(arg[9]);
    polyPRamp = atof(arg[10]); // factor of polymerization/depolimerization with respect to insetion frequency    
    CP_dur[CP_M]=atof(arg[11]);
    polyMRamp = atof(arg[12]); // factor of polymerization/depolimerization with respect to insetion frequency    
    CP_dur[CP_Q]=atof(arg[13]);

    prolRamp = 1.0;
    prolMaxInc  = 0.1;
    
    
    // allocate space to store cell ID's
    nCell = narg-fixargs;
    memory->create(idCell,nCell,"fix/sem:idCell");
    maxIdCell=0;
    for (int i=0; i<nCell; i++){
	idCell[i]=atoi(arg[fixargs+i]);
	maxIdCell=MAX(maxIdCell,idCell[i]);
    }
    memory->create(dtCell,nCell,"fix/sem:dtCell");
    memory->create(idurCell,nCell,"fix/sem:idurCell");
    memory->create(cpCell,nCell,"fix/sem:cpCell");
    memory->create(ccCell,nCell,"fix/sem:ccCell");
    memory->create(tyCell,nCell,"fix/sem:tyCell"); // could be stored globally..
    memory->create(nSceCell,nCell,"fix/sem:nSceCell");
    
    memory->create(polarity,nCell,3,"fix/sem:polarity");
    
    memory->create(bCreateSce,nCell,"fix/sem:createSce");
    memory->create(bSplitCell,nCell,"fix/sem:splitCell");
    memory->create(bPrepareSplit,nCell,"fix/sem:prepareSplit");
    memory->create(bPrepareMig,nCell,"fix/sem:prepareMig");
    memory->create(bMig,nCell,"fix/sem:Mig");
 
    // enable force reneighboring
    force_reneighbor = 1;
    next_reneighbor = -1;
}
/* ----------------------------------------------------------------------
 Destruct
 ---------------------------------------------------------------------- */
FixSemPM::~FixSemPM()
{
    memory->destroy(idCell);
    memory->destroy(dtCell);
    memory->destroy(idurCell);
    memory->destroy(cpCell);
    memory->destroy(ccCell);
    memory->destroy(tyCell);    
    memory->destroy(nSceCell);
    
    memory->destroy(polarity);
    
    memory->destroy(bCreateSce);
    memory->destroy(bSplitCell);
    memory->destroy(bPrepareSplit);
    memory->destroy(bPrepareMig);
    memory->destroy(bMig);
}
/* ----------------------------------------------------------------------
 When to invoce fix during execution pipeline
 ---------------------------------------------------------------------- */
int FixSemPM::setmask()
{
    int mask = 0;
    mask |= PRE_FORCE;
    mask |= END_OF_STEP;
    return mask;
}
/* ----------------------------------------------------------------------
 Init things:
 ---------------------------------------------------------------------- */
void FixSemPM::init()
{
    // check some things
    if (!atom->sem_flag) error->all(FLERR,"Fix adapt requires atom of sem style");
    
    polyPInc=((double)sizeCell)*update->dt*((double)nevery)*polyPRamp/CC_dur[CC_I];
    polyMInc=((double)sizeCell)*update->dt*((double)nevery)*polyMRamp/CP_dur[CP_M];    
#ifdef __PROL__
    preprolInc=update->dt*((double)nevery)*prolRamp*prolMaxInc/(CC_dur[CC_M1]);
    postprolInc=update->dt*((double)nevery)*prolRamp*prolMaxInc/(CC_dur[CC_M2]);
    printf("FixProliferateInitiated: cellSize %d, G0 %f, I %f M1 %f M2 %f polyInc %f preprolInc %f postprolInc %f",
	   sizeCell,CC_dur[CC_G0],CC_dur[CC_I],CC_dur[CC_M1],CC_dur[CC_M2],polyPInc,preprolInc,postprolInc);
    
#endif
    
    int *nSceCellLoc, *tyCellLoc;
    memory->create(nSceCellLoc,nCell,"fix/sem:nSceCellLoc");
    memory->create(tyCellLoc,nCell,"fix/sem:tyCellLoc");
    
    // reset things
    for (int i=0; i<nCell; i++){
	dtCell[i]=0.0;
	ccCell[i]=CC_I; //@@@ maybe set to G0
	//@@@ Sandipan-- change functionality from proliferation to migration .. default is CP_P
	cpCell[i]=CP_P;
	//cpCell[i]=CP_Q;
	//@@@ Sandipan-- change ends
	nSceCellLoc[i]=0;
	tyCellLoc[i]=0;
	polarity[i][0]=1.0;
	polarity[i][1]=0.0;
	polarity[i][2]=0.0;
    }
    int nlocal=atom->nlocal;
    int *cell=atom->cell;
    int *type=atom->type;
    int *flag;
    memory->create(flag,nlocal,"fix/sem:flag");
    
    // init number of sce's per cell
    for(int i=0;i<nlocal;i++){
	flag[i]=0;
	for(int j=0;j<nCell;j++){
	    if(cell[i]==idCell[j]){
		nSceCellLoc[j]++;
		tyCellLoc[j]=type[i];
		break;
	    }
	}
    }    
    MPI_Allreduce(nSceCellLoc,nSceCell,nCell,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(tyCellLoc,tyCell,nCell,MPI_INT,MPI_MAX,world);
    
    //create groups for poly/depoly proliferation
    group->create("depolyP",flag);
    depolyPGroupId=group->find("depolyP");
    depolyPGroupBitmask=group->bitmask[depolyPGroupId];
    group->create("polyP",flag);
    polyPGroupId=group->find("polyP");
    polyPGroupBitmask=group->bitmask[polyPGroupId];
    
    //create groups for poly/depoly migration
    group->create("depolyM",flag);
    depolyMGroupId=group->find("depolyM");
    depolyMGroupBitmask=group->bitmask[depolyMGroupId];
    group->create("polyM",flag);
    polyMGroupId=group->find("polyM");
    polyMGroupBitmask=group->bitmask[polyMGroupId];
    
    //create groups for proliferating cells
    group->create("preprol",flag);
    preprolGroupId=group->find("preprol");
    preprolGroupBitmask=group->bitmask[preprolGroupId];
    group->create("postprol",flag);
    postprolGroupId=group->find("postprol");
    postprolGroupBitmask=group->bitmask[postprolGroupId];
    
    
    memory->destroy(nSceCellLoc);
    memory->destroy(tyCellLoc);
    memory->destroy(flag);
}
/* ---------------------------------------------------------------------
 - advance cell clock
 - create element
 - advance cycle if necessary
 ---------------------------------------------------------------------- */
void FixSemPM::final_integrate(){
    
}

/* ----------------------------------------------------------------------
 What to do after step
 ---------------------------------------------------------------------- */
void FixSemPM::end_of_step(){
    if (nevery == 0) return;
    if (update->ntimestep % nevery) return;
    advanceCell(); 
}


/* ----------------------------------------------------------------------
 What to do before forces are calculated
 ---------------------------------------------------------------------- */
void FixSemPM::pre_force(int){
    
    //reset cell dens
    double * dens = atom->dens;
    int nlocal = atom->nlocal;
    for(int i=0;i<nlocal;i++) dens[i]=0.0;
    
    
}
/* ----------------------------------------------------------------------
 What to do after run
 ---------------------------------------------------------------------- */
void FixSemPM::post_run(){
    
}
/* ----------------------------------------------------------------------
 for all cells
 - advance cell clock
 - ?create element
 - ?switch cell cycle phase
 
 ---------------------------------------------------------------------- */
void FixSemPM::advanceCell(){
#ifdef __SEM_DEBUG__
    printf(BIGINT_FORMAT " advanceCell\n",update->ntimestep);
#endif
    int me;
    MPI_Comm_rank(world,&me);
    
    // init params
    nCreateSce=nSplitCell=nPrepareSplit=nPrepareMig=nMig=0;
    
    int nProl=0;
    
    //++++++++++++++++++++++++++++++++++
    //	- run through all cells
    //	    - advance cell clock
    //	    - check cell phenotype
    //		- register polymerization events
    //
    //  - do the polymerization
    //  - do the splitting
    //++++++++++++++++++++++++++++++++++
    for (int i=0;i<nCell;i++){
	bCreateSce[i]=false;
	bPrepareSplit[i]=false;
	bSplitCell[i]=false;
	bPrepareMig[i]=false;
	bMig[i]=false;
	dtCell[i]+=update->dt*(double)nevery;
	switch (cpCell[i]){ // CELL PHENOTYPE CHECK
	    case CP_Q: // QUIESCENT CELL
		// do nothing...
		//@@@ for the moment, do nothing for 
		//@@@ then switch to migration;
		if(dtCell[i]>CP_dur[CP_Q]){
		    dtCell[i]-=CP_dur[CP_Q];
		    bPrepareMig[i]=true;
		    nPrepareMig++;
		    double fluc = 0.0;
		    idurCell[i]=CP_dur[CP_M]/(double)sizeCell+fluc;
		    cpCell[i]=CP_M;
		    printf(BIGINT_FORMAT ": cell %d switch to CP_M: dtCell=%f, nSce=%d\n",update->ntimestep,idCell[i],dtCell[i],nSceCell[i]);
		}
		break;
	    case CP_M: // MIGRATING CELL
		// coordinated polymerization/depolymerization
		// at defined intervals
		if(dtCell[i]>idurCell[i]){
		    // register for migration polymerization
		    dtCell[i]-=idurCell[i];
		    bMig[i]=true;
		    nMig++;
		    double fluc=0.0;
		    idurCell[i]=CP_dur[CP_M]/(double)sizeCell+fluc;
		}
		break;
	    case CP_P: // PROLIFERATING CELL
		// run cell cycle dependent proliferation program
		switch (ccCell[i]){ // CELL CYCLE CHECK
		    case CC_G0: // RESTING PHASE
			if(dtCell[i]>CC_dur[CC_G0]){
			    // advance in cell cycle phase
			    dtCell[i]-=CC_dur[CC_G0];
			    double fluc=0.0;
			    idurCell[i]=CC_dur[CC_I]/(double)sizeCell+fluc;
			    ccCell[i]=CC_I;
			    printf(BIGINT_FORMAT ": cell %d switch to CC_I: dtCell=%f, nSce=%d\n",update->ntimestep,idCell[i],dtCell[i],nSceCell[i]);
			}else{
			    break;
			}
		    case CC_I: // GROWTH PHASE
			if(dtCell[i]>idurCell[i]){
			    dtCell[i]-=idurCell[i];
			    // check if advance in cell cycle phase
			    
			    //@@@ FIX: grow and stop
			    if(nSceCell[i]==sizeCell){
				// jump quiescent
				cpCell[i]=CP_Q;
				break;
			    }
			    //@@@ END FIX
			    if(nSceCell[i]<sizeCell*2){
				// create sce
				bCreateSce[i]=true;
				nCreateSce++;
				double fluc=0.0;
				idurCell[i]=CC_dur[CC_I]/(double)sizeCell+fluc;			
#ifdef __SEM_DEBUG__
				printf(BIGINT_FORMAT ": cell %d creating SCE: dtCell=%f, nSce=%d\n",update->ntimestep,idCell[i],dtCell[i],nSceCell[i]);
#endif
				break;
			    }else{
				bPrepareSplit[i]=true;
				nPrepareSplit++;
				ccCell[i]=CC_M1;
				printf(BIGINT_FORMAT ": cell %d switch to CC_M1: dtCell=%f, nSce=%d\n",update->ntimestep,idCell[i],dtCell[i],nSceCell[i]);
			    }
			}else{
			    break;
			}
		    case CC_M1: // SEPARATION PHASE
			// wait a bit before dividing
			if(dtCell[i]>=CC_dur[CC_M1]){
			    dtCell[i]-=CC_dur[CC_M1];
			    ccCell[i]=CC_M2;
			    bSplitCell[i]=true;
			    nSplitCell++;
			    printf(BIGINT_FORMAT ": cell %i switch to CC_M2: dtCell=%f, nSce=%d\n",update->ntimestep,idCell[i],dtCell[i],nSceCell[i]);
			}else{
#ifdef __PROL__
			    nProl++;
#endif
			    break;
			}
		    case CC_M2: // SPLITTING PHASE (post splitting
			// wait a bit after dividing
			if(dtCell[i]>=CC_dur[CC_M2]){
			    dtCell[i]-=CC_dur[CC_M2];
			    ccCell[i]=CC_G0;
			    printf(BIGINT_FORMAT ": cell %i switch to CC_G0: dtCell=%f, nSce=%d\n",update->ntimestep,idCell[i],dtCell[i],nSceCell[i]);
			    break;
			}else{
			    break;
#ifdef __PROL__
			    nProl++;
#endif		    
			}
		
		    default:
			error->all(FLERR,"Illegal Cell Cycle Phase");
		}
		break;
	    default:
		error->all(FLERR,"Illegal Cell Phenotype");
	}
    }
    
    //++++++++++++++++++++++++++++++++++
    // DO THINGS
    //	- POLYMERIZATION
    //	- SPLITTING
    //  - MIGRATION
    //++++++++++++++++++++++++++++++++++
    if(nCreateSce>0)createSce();
    if(nMig>0)migrate();
	
    if(nPolymerizeSce>0 || nProl>0)polymerizeSce();
#ifdef __PROL__
    if(nPrepareSplit>0) prepareSplit();
#endif
    if(nSplitCell>0) splitCell();
    if(nPrepareMig>0) prepareMig();
    

#ifdef __SEM_DEBUG__
    printf(BIGINT_FORMAT " doneAdvanceCell\n",update->ntimestep);
#endif
}

/* ----------------------------------------------------------------------
    choose call one of the createPolymerizeSce methods
 ---------------------------------------------------------------------- */
void FixSemPM::createSce(){
  // createSceCenter();
    createSceRandom();
}

/* ----------------------------------------------------------------------
    - create SCE for tagged cells at the center of the cell
    - add SCE to polymerize group
 ---------------------------------------------------------------------- */
void FixSemPM::createSceCenter()
{
#ifdef __SEM_DEBUG__
    printf(" - Creating %d SCEs Center\n",nCreateSce);
#endif
    int nlocal;
    double **selCenter, **selCenterLoc;
    int * selNSce,*selId,*selInd;
    
    // init things
    int nSel=nCreateSce;
    bool* bSel=bCreateSce;
    memory->create(selCenter,nSel,3,"fix/sem:selCenter");
    memory->create(selCenterLoc,nSel,3,"fix/sem:selCenterLoc");
    memory->create(selNSce,nSel,"fix/sem:selNSce");
    memory->create(selId,nSel,"fix/sem:selId");
    memory->create(selInd,nSel,"fix/sem:selInd");
    
    int cur=0;
    for(int c=0;c<nCell;c++){
	if(bSel[c]){
	    selCenterLoc[cur][0]=0.0;
	    selCenterLoc[cur][1]=0.0;
	    selCenterLoc[cur][2]=0.0;
	    selNSce[cur]=nSceCell[c];
	    selId[cur]=idCell[c];
	    selInd[cur]=c;
	    cur++;
	}
    }
    
    // find centers of growing cells    
    // accumulate locations localy
    int me;
    MPI_Comm_rank(world,&me);
    int cnt=0;
    
    nlocal = atom->nlocal;
    int *cell = atom->cell;
    int *mask  = atom->mask;
    double **x = atom->x;
    double *p = atom->p;
    
    for (int i=0;i<nlocal;i++){	
	for(int j=0;j<nSel;j++){
	    if(cell[i]==selId[j]){ // @@@??? maybe just do for all
		
		selCenterLoc[j][0]+=x[i][0];
		selCenterLoc[j][1]+=x[i][1];
		selCenterLoc[j][2]+=x[i][2];
		cnt++;
	    }
	}
    }
    //sum up locations globally
    
    //    printf("PROC: %d cnt %d\n",me,cnt);
    for(int i=0;i<nSel;i++){
	//	printf("proc: %d cell: %d | centerLoc: %f %f %f | selNSce %d \n",me,selId[i],selCenterLoc[i][0],selCenterLoc[i][1],selCenterLoc[i][2], selNSce[i]);
    }
    MPI_Allreduce(&selCenterLoc[0][0],&selCenter[0][0],3*nSel,MPI_DOUBLE,MPI_SUM,world);

    
    for(int i=0;i<nSel;i++){
	//	printf("proc: %d cell: %d | center: %f %f %f\n",me,selId[i],selCenter[i][0],selCenter[i][1],selCenter[i][2]);
	selCenter[i][0]/=(double)selNSce[i];
	selCenter[i][1]/=(double)selNSce[i];
	selCenter[i][2]/=(double)selNSce[i];
    }
    
    
    // create atoms
    // - check if in my box
    // - set tag
    // - invoce set_arrays() fir fixes that need initialization of new atoms
    bigint natoms_previous = atom->natoms;
    int nlocal_previous = atom->nlocal;
    double lambda[3],*coord;
    double *sublo,*subhi;
    
    if (domain->triclinic == 0){
	sublo = domain->sublo;
	subhi = domain->subhi;
    }else{
	printf("SEM-WARNING:triclinic not supported but not prohibited.. maybe it works...");
	sublo = domain->sublo_lamda;
	subhi = domain->subhi_lamda;
    }
    
    //    printf("proc %d, sublo %f %f %f, subhi %f %f %f\n",me,sublo[0],sublo[1],sublo[2],subhi[0],subhi[1],subhi[2]);
    for(int i=0;i<nSel;i++){
	
	if (domain->triclinic == 0){
	    coord = selCenter[i];
	}else{
	    printf("SEM-WARNING:triclinic not supported but not prohibited.. maybe it works...");
	    domain->x2lamda(selCenter[i],lambda);
	    coord = lambda;
	}
	// create atom if in my subbox
	if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
	    coord[1] >= sublo[1] && coord[1] < subhi[1] &&
	    coord[2] >= sublo[2] && coord[2] < subhi[2]){
#ifdef __SEM_DEBUG__
	    printf("--" BIGINT_FORMAT " | proc: %d cell: %d  CREATING IT center(%f,%f,%f)\n",update->ntimestep,me,selId[i],selCenter[i][0],selCenter[i][1],selCenter[i][2]);
	    
#endif
	    
	    atom->avec->create_atom(tyCell[selInd[i]],coord);
	    atom->cell[atom->nlocal-1]=selId[i];
	    atom->p[atom->nlocal-1]=0.0;
	    atom->mask[atom->nlocal-1] |= polyPGroupBitmask;
	    
#ifdef __SEM_DEBUG__
	    printf("id %d, center(%f,%f,%f), type %d, cell %d, p %f,\n",atom->tag[atom->nlocal-1],atom->x[atom->nlocal-1][0],atom->x[atom->nlocal-1][1],atom->x[atom->nlocal-1][2],atom->type[atom->nlocal-1],atom->cell[atom->nlocal-1],atom->p[atom->nlocal-1]);
#endif
	}
	nSceCell[selInd[i]]++;
	nPolymerizeSce++;
    }
    // invoke set_arrays() for fixes that need initialization of new atoms
    nlocal = atom->nlocal;
    for (int m = 0; m < modify->nfix; m++) {
	Fix *fix = modify->fix[m];
	if (fix->create_attribute)
	    for (int i = nlocal_previous; i < nlocal; i++)
		fix->set_arrays(i);
    }
    
    // new total # of atoms
    bigint rlocal = atom->nlocal;
    MPI_Allreduce(&rlocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
    if (atom->natoms < 0 || atom->natoms > MAXBIGINT)
	error->all(FLERR,"Too many total atoms");
    
    
    // @@@SEM not sure about this?
    // set tag # of new particles beyond all previous atoms
    // if global map exists, reset it now instead of waiting for comm
    // since deleting atoms messes up ghosts
    if (atom->tag_enable) {
      atom->tag_extend();
      if (atom->map_style) {
	atom->nghost = 0;
	atom->map_init();
	atom->map_set();
      }
    }
    
    // never invoced
    if (atom->molecular) {
	int **nspecial = atom->nspecial;
	for (int i = nlocal_previous; i < atom->nlocal; i++) {
	    nspecial[i][0] = 0;
	    nspecial[i][1] = 0;
	    nspecial[i][2] = 0;
	}
    }
    
    //@@@SEM: force rebuild of neighbour List next time step:
    //needs neighbour must_check enabled!!!
    next_reneighbor=update->ntimestep+1;
    
    
    //    printf("proc: %d natoms: %d, nlocal %d \n",me, atom->natoms,atom->nlocal);
    // free memory
    memory->destroy(selCenter);
    memory->destroy(selCenterLoc);
    memory->destroy(selNSce);
    memory->destroy(selId);
    memory->destroy(selInd);
    
}
    
/* ----------------------------------------------------------------------
 - create SCE for tagged cells by duplicating random element
    - select a random element from nSceCell
    - rank 0 handles the shizzl
	- count nSceCellLoc for growing cells
	- send to rank 0
	- 
 - add SCE to polymerize group
---------------------------------------------------------------------- */
void FixSemPM::createSceRandom()
{
#ifdef __SEM_DEBUG__
    printf(" - Creating %d SCEs Random\n",nCreateSce);
#endif
    
    //INIT THINGS
    int me;
    MPI_Comm_rank(world,&me);
    int nlocal;
    int * selNSceLoc,*selId,*selInd;
    int *selElemProc;
    
    // only allocated by rank 0
    int * selNSceAll, * selNSce;
	    
    // init things
    int nSel=nCreateSce;
    bool* bSel=bCreateSce;
    
    memory->create(selNSceLoc,nSel,"fix/sem:selNSceLoc");
    memory->create(selId,nSel,"fix/sem:selId");
    memory->create(selInd,nSel,"fix/sem:selInd");
    memory->create(selElemProc,2*nSel,"fix/sem:selNElemProc");
    
    // init to 0
    int cur=0;
    for(int c=0;c<nCell;c++){
	if(bSel[c]){
	    selNSceLoc[cur]=0.0;
	    selId[cur]=idCell[c];
	    selInd[cur]=c;
	    cur++;
	    
	    // doin it here.. should be done after or when created...
	    nSceCell[c]++;
	}
    }
    
    // count number of SCEs on this proc    
    nlocal = atom->nlocal;
    int *cell = atom->cell;
    double **x = atom->x;
    int *mask  = atom->mask;
    double *p = atom->p;
    int *cellpart = atom->cellpart;
    
    for (int i=0;i<nlocal;i++){
	for(int j=0;j<nSel;j++){
	    if(cell[i]==selId[j]&&cellpart[i]!=2){ // @@@maybe exclude membrane elements
		selNSceLoc[j]++;
	    }
	}
    }
    
#ifdef __SEM_DEBUG__
    for(int i=0;i<nSel;i++){
	printf("id: %d, selNSceLoc[%d]: %d\n",me,i,selNSceLoc[i]);
    }
    printf("\n");
#endif

    nPolymerizeSce+=nSel;
    
    //Gather all information on master node
    // follow: http://www.open-mpi.org/doc/v1.5/man3/MPI_Gather.3.php
    int gsize;
    if( me == 0){
	MPI_Comm_size(world,&gsize);
	// allocate arrays needed on master
        memory->create(selNSceAll,gsize*nSel,"fix/sem:selNSceAll");
        memory->create(selNSce,nSel,"fix/sem:selNSce");
    }
    MPI_Gather(selNSceLoc,nSel,MPI_INT,selNSceAll,nSel, MPI_INT, 0,world);


    // on rank 0, decide who to duplicate element
    if( me == 0){
#ifdef __SEM_DEBUG__
	for(int i=0;i<nSel*gsize;i++){
	    printf("id: %d, selNsceAll[%d]: %d\n",me,i,selNSceAll[i]);
	}
	printf("\n");
#endif	
	// count elements per cell
	// pick random element
	for(int i=0;i<nSel;i++){
	    selNSce[i]=0;
	    for(int j=0;j<gsize;j++){
		selNSce[i]+=selNSceAll[j*nSel+i];
	    }
	    selElemProc[i*2]=rand()%selNSce[i];
	    selElemProc[i*2+1]=0;
	}
	
#ifdef __SEM_DEBUG__
	for(int i=0;i<nSel*2;i++){
	    printf("id: %d, selElemProc[%d]: %d\n",me,i,selElemProc[i]);
	}
	printf("\n");
#endif
	
	//map picked elements to element on processor
	for (int i=0;i<nSel;i++){
	    for(int j=0;j<gsize;j++){
		if(selElemProc[2*i]>=selNSceAll[j*nSel+i]){
		    selElemProc[2*i]-=selNSceAll[j*nSel+i];
		}else{
		    selElemProc[2*i+1]=j;
		    break;
		}
	    }
	}
#ifdef __SEM_DEBUG__
	for(int i=0;i<nSel*2;i++){
	    printf("id: %d, selElemProc[%d]: %d\n",me,i,selElemProc[i]);
	}
	printf("\n");
#endif
	
	//free master memory
	memory->destroy(selNSceAll);
	memory->destroy(selNSce);
    }
    // broadcast to all
    MPI_Bcast(selElemProc,nSel*2,MPI_INT,0,world);
	      
    
#ifdef __SEM_DEBUG__
    printf("AFTERMPI\n");
    for(int i=0;i<nSel*2;i++){
	printf("id: %d, selElemProc[%d]: %d\n",me,i,selElemProc[i]);
    }
    printf("\n");
#endif
    
    // count how many I am responsible for
    int nCreateLoc=0;
    for(int i=0;i<nSel;i++){
	if(selElemProc[2*i+1]==me) nCreateLoc++;
    }
    
    
    // some things for after sce's are created
    bigint natoms_previous = atom->natoms;
    int nlocal_previous = atom->nlocal;

    
    if(nCreateLoc>0){
#ifdef __SEM_DEBUG__
	    printf("id %d, creating %d SCEs\n",me, nCreateLoc);
#endif
	// allocate stuff
	int * dupInd;
	int * dupId;
	int * dupElem;
        double ** dupCoord;
        
        memory->create(dupInd,nCreateLoc,"fix/sem:dupInd");
        memory->create(dupId,nCreateLoc,"fix/sem:dupId");
        memory->create(dupElem,nCreateLoc,"fix/sem:dupElem");
        memory->create(dupCoord,nCreateLoc,3,"fix/sem:dupCoord");
	
	// register which ones i am responsible for
	int cnt=0;
	for(int i=0;i<nSel;i++){
	    if(selElemProc[2*i+1]==me){
		dupInd[cnt]=selInd[i];
		dupId[cnt]=selId[i];
		dupElem[cnt]=selElemProc[2*i];
		cnt++;
	    }
	}
	
	// run through elements and store duplication location
	cnt=0;
	
	nlocal = atom->nlocal;
	for(int i=0;i<nlocal;i++){
	    for(int j=0;j<nCreateLoc;j++){
		if(cell[i]==dupId[j]&&cellpart[i]!=2){
		    if(dupElem[j]==0){
			// duplicate this one
			dupCoord[j][0]=x[i][0];
			dupCoord[j][1]=x[i][1];
			dupCoord[j][2]=x[i][2];
			dupElem[j]--;
			cnt++;
			if(cnt==nCreateLoc) break;
		    }else{
			dupElem[j]--;
		    }
		}
		
	    }
	    if(cnt==nCreateLoc) break;
	}
	
	// create elements
	double lambda[3],*coord;
	for(int i=0;i<nCreateLoc;i++){
	    // transform coords
	    if (domain->triclinic == 0){
		coord = dupCoord[i];
	    }else{
		printf("SEM-WARNING:triclinic not supported but not prohibited.. maybe it works...");
		domain->x2lamda(dupCoord[i],lambda);
		coord = lambda;
	    }
	    
	    // create element
	    atom->avec->create_atom(tyCell[dupInd[i]],coord);
	    atom->cell[atom->nlocal-1]=idCell[dupInd[i]];//dupId[i];
	    atom->p[atom->nlocal-1]=0.0;
	    atom->mask[atom->nlocal-1] |= polyPGroupBitmask;
#ifdef __SEM_DEBUG__
	    printf("CREATED IT, cell id = %d %d %d \n",atom->cell[atom->nlocal-1],idCell[dupInd[i]],dupId[i]);
#endif

	}
	
	//dealocate local arrays
	memory->destroy(dupInd);
	memory->destroy(dupId);
	memory->destroy(dupElem);
	memory->destroy(dupCoord);
    }
	    
    // invoke set_arrays() for fixes that need initialization of new atoms    
    nlocal = atom->nlocal;
    
#ifdef __SEM_DEBUG__
    printf("id %d: nlocal before create: %d\n",me, nlocal_previous);
    printf("id %d: nlocal after create: %d\n",me, nlocal);
#endif

    for (int m = 0; m < modify->nfix; m++) {
	Fix *fix = modify->fix[m];
	if (fix->create_attribute){
	    for (int i = nlocal_previous; i < nlocal; i++){
		printf("EVER CALLED?");
		fix->set_arrays(i);
	    }
	}
    }
    
    // new total # of atoms
    bigint rlocal = atom->nlocal;
    MPI_Allreduce(&rlocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
    if (atom->natoms < 0 || atom->natoms > MAXBIGINT)
	error->all(FLERR,"Too many total atoms");
    
	    
    // @@@SEM not sure about this?
    // set tag # of new particles beyond all previous atoms
    // if global map exists, reset it now instead of waiting for comm
    // since deleting atoms messes up ghosts
    if (atom->tag_enable) {
      atom->tag_extend();
      if (atom->map_style) {
	atom->nghost = 0;
	atom->map_init();
	atom->map_set();
      }
    }
    
    // never invoced
    if (atom->molecular) {
	int **nspecial = atom->nspecial;
	for (int i = nlocal_previous; i < atom->nlocal; i++) {
	    nspecial[i][0] = 0;
	    nspecial[i][1] = 0;
	    nspecial[i][2] = 0;
	}
    }
    
    //@@@SEM: force rebuild of neighbour List next time step:
    //needs neighbour must_check enabled!!!
    next_reneighbor=update->ntimestep+1;
    
	    
    //    printf("proc: %d natoms: %d, nlocal %d \n",me, atom->natoms,atom->nlocal);
    // free memory    
    memory->destroy(selNSceLoc);
    memory->destroy(selId);
    memory->destroy(selInd);
    memory->destroy(selElemProc);
    MPI_Barrier(world); //@@@ remove!!
}
	
/* ----------------------------------------------------------------------
 polymerize in group
 ---------------------------------------------------------------------- */
void FixSemPM::polymerizeSce(){
#ifdef __SEM_DEBUG__
    printf(" - Polymerizing %d SCEs\n",nPolymerizeSce);
#endif
    int nPolymerizeSceLoc=0;
    int nlocal = atom->nlocal;
    int *cell = atom->cell;
    int *mask  = atom->mask;
    double **x = atom->x;
    double *p = atom->p;
    
    int me;
    MPI_Comm_rank(world,&me);
    
    int nKilled,nKilledLoc;
    nKilledLoc=0;
    
    // go through backwards to delete elements
    for (int i=nlocal-1;i>=0;i--){
	//	printf("--" BIGINT_FORMAT " | proc: %d, particle %d, location: %f %f %f\n",update->ntimestep, me, atom->tag[i],x[i][0],x[i][1],x[i][2]);
	// POLYERIZATION PROLIFERATION
	if(mask[i] & polyPGroupBitmask){
	    nPolymerizeSceLoc++;
	    p[i]+=polyPInc;
	    if(p[i]>1.0){
		p[i]=1.0;
		mask[i]^=polyPGroupBitmask;
		nPolymerizeSceLoc--;
	    }
	}
	// DEPOLYMERIZATION PROLIFERATION
	else if(mask[i] & depolyPGroupBitmask){
	    if(p[i]-=polyPInc<0.0){		
		//p[i]=0;
		//mask[i]^=depolyPGroupBitmask;
		nPolymerizeSceLoc--;
		atom->avec->copy(atom->nlocal-1,i,0);
		atom->nlocal--;
		nKilledLoc++;
	    }
	}
	// POLYMERIZATION MIGRATION
	else if(mask[i] & polyMGroupBitmask){
	    nPolymerizeSceLoc++;
	    p[i]+=polyMInc;
	    if(p[i]>1.0){
		p[i]=1.0;
		mask[i]^=polyMGroupBitmask;
		nPolymerizeSceLoc--;
	    }
	}
	// DEPOLYMERIZATION MIGRATION
	else if(mask[i] & depolyMGroupBitmask){
	    p[i]-=polyMInc;
	    if(p[i]<0.0){
		//p[i]=0;
		//mask[i]^=depolyMGroupBitmask;
		nPolymerizeSceLoc--;
		atom->avec->copy(atom->nlocal-1,i,0);
		atom->nlocal--;
		nKilledLoc++;
	    }
	}
#ifdef __PROL__
	// BLOW/SHRINK 4 PROLIFERATION
	// 
	else{
	    if(mask[i] & preprolGroupBitmask){
		p[i]+=preprolInc;
	    }else if(mask[i] & postprolGroupBitmask){
		if(p[i]-=postprolInc<1.0){
		    p[i]=1.0;
		    mask[i]^=postprolGroupBitmask;
		}
	    }
	}
#endif
    }
    
    MPI_Allreduce(&nKilledLoc,&nKilled,1,MPI_INT,MPI_SUM,world);
    if(nKilled>0){
#ifdef __SEM_DEBUG__
	printf("killed: %d\n",nKilled);
#endif
	atom->natoms -= nKilled;
	//@@@SEM: force rebuild of neighbour List next time step:
	//needs neighbour must_check enabled!!!
	next_reneighbor=update->ntimestep+1;
    }
    // update atom->natoms;
    MPI_Allreduce(&nPolymerizeSceLoc,&nPolymerizeSce,1,MPI_INT,MPI_SUM,world);
}


/* ----------------------------------------------------------------------
 prepare splitting
 - tag elements of cell that is about to split with preprolGroupBitmap
 ---------------------------------------------------------------------- */
void FixSemPM::prepareSplit(){
#ifdef __SEM_DEBUG__
    printf(" - prepareSplitting %d Cells of %d Cells in total\n",nPrepareSplit,nCell);
#endif
    int *selId;
    int nSel;
    bool *bSel;
    
    nSel=nPrepareSplit;
    bSel=bPrepareSplit;
    memory->create(selId,nSel,"fix/sem:selId");
    
    int cur=0;
    for(int i=0;i<nCell;i++){
	if(bSel[i]){
	    selId[cur]=idCell[i];
	    cur++;
	}	
    }
    
    int nlocal=atom->nlocal;
    int* cell=atom->cell;
    int* mask=atom->mask;
    for(int i=0;i<nlocal;i++){
	for(int j=0;j<nSel;j++){
	    if(cell[i]==selId[j])mask[i]|=preprolGroupBitmask;
	}
    }	
    
    memory->destroy(selId);
}

/* ----------------------------------------------------------------------
 split tagged cells
 - find center [x0]
 - find R [sum(x-x0)*(x-x0)]
 ---------------------------------------------------------------------- */
void FixSemPM::splitCell(){
#ifdef __SEM_DEBUG__
    printf(" - Splitting %d Cells of %d Cells in total\n",nSplitCell,nCell);
    int me;
    MPI_Comm_rank(world,&me);
#endif
    int nlocal;
    double **selCenter, **selCenterLoc;
    double (* selEv)[3][3];
    double (* selR)[3][3],(* selRLoc)[3][3];
    double (* selEw)[3];
    int * selNSce, *selNSceLoc,*selId,*selInd,*selIndE;
    
    // init things
    int nSel=nSplitCell;
    
    memory->create(selCenter,nSel,3,"fix/sem:selCenter");
    memory->create(selCenterLoc,nSel,3,"fix/sem:selCenterLoc");
    memory->create(selNSce,nSel,"fix/sem:selNSce");
    memory->create(selNSceLoc,nSel,"fix/sem:selNSceLoc");
    memory->create(selId,nSel,"fix/sem:selId");
    memory->create(selInd,nSel,"fix/sem:selInd"); 
    
    memory->create(selR,nSel,"fix/sem:selR");
    memory->create(selRLoc,nSel,"fix/sem:selRLoc");
    memory->create(selEv,nSel,"fix/sem:selEv");
    memory->create(selEw,nSel,"fix/sem:selEw");
    memory->create(selIndE,nSel,"fix/sem:selIndE");
    
    //+++++++++++++++++
    // extend arrays to support new cells
    //+++++++++++++++++
#ifdef __SEM_DEBUG__
    printf("before array extend \n");
    for(int i=0;i<nCell;i++){
	printf("Cell %d: idCell %d, dtCell %f, idurCell %f, ccCell %d, tyCell %d, nSceCell %d, bCreateSce %d, bSplitCell %d\n", 
	       i, idCell[i], dtCell[i], idurCell[i], ccCell[i], tyCell[i], nSceCell[i], bCreateSce[i], bSplitCell[i]);
	if(idCell[i]>nCell || idCell[i]<0)	error->all(FLERR,"Illegal cell id 1");
    }
#endif
    int nCellNew=nCell+nSel;
    memory->grow(idCell,nCellNew,"fix/sem:idCell");
    memory->grow(dtCell,nCellNew,"fix/sem:dtCell");
    memory->grow(idurCell,nCellNew,"fix/sem:idurCell");
    memory->grow(ccCell,nCellNew,"fix/sem:ccCell");
    memory->grow(cpCell,nCellNew,"fix/sem:cpCell");
    memory->grow(tyCell,nCellNew,"fix/sem:tyCell"); // could be stored globally..
    memory->grow(nSceCell,nCellNew,"fix/sem:nSceCell");
    
    memory->grow(polarity,nCellNew,3,"fix/sem:polarity");
    
    memory->grow(bCreateSce,nCellNew,"fix/sem:bCreateSce");
    memory->grow(bSplitCell,nCellNew,"fix/sem:bSplitCell");
    memory->grow(bPrepareSplit,nCellNew,"fix/sem:bPrepareSplit");
    memory->grow(bPrepareMig,nCellNew,"fix/sem:bPrepareMig");
    memory->grow(bMig,nCellNew,"fix/sem:bMig");
    
    bool* bSel=bSplitCell;
    
    
    //+++++++++++++++++
    // find Cell Center
    //+++++++++++++++++
    
    int cur=0;
    for(int c=0;c<nCell;c++){
	if(bSel[c]){
	    selCenterLoc[cur][0]=0.0;
	    selCenterLoc[cur][1]=0.0;
	    selCenterLoc[cur][2]=0.0;
	    selNSce[cur]=nSceCell[c];
	    selNSceLoc[cur]=0;
	    selId[cur]=idCell[c];
	    selInd[cur]=c;
	    for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
		    selR[cur][i][j]=0.0;
	    // init extended arrays
	    idCell[nCell+cur]=++maxIdCell;
	    dtCell[nCell+cur]=dtCell[c];
	    idurCell[nCell+cur]=idurCell[c];
	    ccCell[nCell+cur]=ccCell[c];
	    cpCell[nCell+cur]=cpCell[c];
	    tyCell[nCell+cur]=tyCell[c];
	    nSceCell[nCell+cur]=0;
	    
	    polarity[nCell+cur][0]=1.0;
	    polarity[nCell+cur][1]=0.0;
	    polarity[nCell+cur][2]=0.0;
	    
	    bCreateSce[nCell+cur]=false;
	    bSplitCell[nCell+cur]=false;
	    bPrepareSplit[nCell+cur]=false;
	    bPrepareMig[nCell+cur]=false;
	    bMig[nCell+cur]=false;
	    cur++;
	}
    }
    
#ifdef __SEM_DEBUG__
    printf("array extend \n");
    for(int i=0;i<nCellNew;i++){
	printf("Cell %d: idCell %d, dtCell %f, idurCell %f, ccCell %d, tyCell %d, nSceCell %d, bCreateSce %d, bSplitCell %d\n", 
	       i, idCell[i], dtCell[i], idurCell[i], ccCell[i], tyCell[i], nSceCell[i], bCreateSce[i], bSplitCell[i]);
	if(idCell[i]>nCellNew || idCell[i]<0)	error->all(FLERR,"Illegal cell id");
    }
#endif
    // find centers of growing cells    
    // accumulate locations localy
    nlocal = atom->nlocal;
    int *cell = atom->cell;
    int *mask  = atom->mask;
    double **x = atom->x;
    double *p = atom->p;
    for (int i=0;i<nlocal;i++){
	for(int j=0;j<nSel;j++){
	    if(cell[i]==selId[j]){ // @@@??? maybe just do it for all...
		selCenterLoc[j][0]+=x[i][0]/(double)selNSce[j];
		selCenterLoc[j][1]+=x[i][1]/(double)selNSce[j];
		selCenterLoc[j][2]+=x[i][2]/(double)selNSce[j];
	    }
	}
    }
#ifdef __SEM_DEBUG__
    printf("comm center... \n");
#endif
    
    
    // sum up locations globally
    MPI_Allreduce(&selCenterLoc[0][0],&selCenter[0][0],3*nSel,MPI_DOUBLE,MPI_SUM,world);
    
#ifdef __SEM_DEBUG__
    printf("accumulate R \n");
#endif
    
    //+++++++++++++++++
    // accumulate R
    //+++++++++++++++++
    for (int i=0;i<nlocal;i++){
	for(int j=0;j<nSel;j++){
	    if(cell[i]==selId[j]){ // @@@??? maybe just do it for all...
		selRLoc[j][0][0]+=(x[i][0]-selCenter[j][0])*(x[i][0]-selCenter[j][0])/(double)selNSce[j];
		selRLoc[j][0][1]+=(x[i][0]-selCenter[j][0])*(x[i][1]-selCenter[j][1])/(double)selNSce[j];
		selRLoc[j][0][2]+=(x[i][0]-selCenter[j][0])*(x[i][2]-selCenter[j][2])/(double)selNSce[j];
		selRLoc[j][1][1]+=(x[i][1]-selCenter[j][1])*(x[i][1]-selCenter[j][1])/(double)selNSce[j];
		selRLoc[j][1][2]+=(x[i][1]-selCenter[j][1])*(x[i][2]-selCenter[j][2])/(double)selNSce[j];
		selRLoc[j][2][2]+=(x[i][2]-selCenter[j][2])*(x[i][2]-selCenter[j][2])/(double)selNSce[j];
	    }
	}
    }
    for(int j=0;j<nSel;j++){
	selRLoc[j][1][0]=selRLoc[j][0][1];
	selRLoc[j][2][0]=selRLoc[j][0][2];
	selRLoc[j][2][1]=selRLoc[j][1][2];
    }
#ifdef __SEM_DEBUG__
    printf("comm R... \n");
#endif
    
    MPI_Allreduce(&selRLoc[0],&selR[0],3*3*nSel,MPI_DOUBLE,MPI_SUM,world);
    
    //+++++++++++++++++
    // calculate Ev/Ew
    //+++++++++++++++++
    for(int j=0;j<nSel;j++){
	eigen_decomposition(selR[j],selEv[j],selEw[j]);
	// largest eigenvalue
	
	// DIVIDE perp to Max EW
	int ind=(fabs(selEw[j][1])>fabs(selEw[j][0]))?1:0;
	selIndE[j]=(fabs(selEw[j][2])>fabs(selEw[j][ind]))?2:ind;
#ifdef __SEM_DEBUG__
	printf("Ew: %f, %f, %f\n",selEw[j][0],selEw[j][1],selEw[j][2]);
	printf("index: %d\n",selIndE[j]);
#endif
	
    }
    
    //+++++++++++++++++
    // reassign sces to new Cell
    //+++++++++++++++++
    for (int i=0;i<nlocal;i++){
	for(int j=0;j<nSel;j++){
	    if(cell[i]==selId[j]){
		// calculate the projection
		double proj=(x[i][0]-selCenter[j][0])*selEv[j][selIndE[j]][0]+
		(x[i][1]-selCenter[j][1])*selEv[j][selIndE[j]][1]+
		(x[i][2]-selCenter[j][2])*selEv[j][selIndE[j]][2];
		if(proj>0){
		    cell[i]=idCell[nCell+j];	// assign particle to new cell
		}else {
		    selNSceLoc[j]++;		// dec sce count
		}
#ifdef __PROL__
		//also enter post proliferation program
		mask[i]^=preprolGroupBitmask;
		mask[i]|=postprolGroupBitmask;
#endif
		
		
	    }
	}
    }
    MPI_Allreduce(&selNSceLoc[0],&selNSce[0],nSel,MPI_INT,MPI_SUM,world);
    for(int i=0;i<nSel;i++){
	int ind = selInd[i];
	nSceCell[nCell+i]=nSceCell[ind]-selNSce[i];
	nSceCell[ind]=selNSce[i];
    }
    
    
    nCell=nCellNew;
#ifdef __SEM_DEBUG__
    printf("Cell Splitted\n");
    for(int i=0;i<nCellNew;i++){
	printf("Cell %d: idCell %d, dtCell %f, idurCell %f, ccCell %d, tyCell %d, nSceCell %d, bCreateSce %d, bSplitCell %d\n", 
	       i, idCell[i], dtCell[i], idurCell[i], ccCell[i], tyCell[i], nSceCell[i], bCreateSce[i], bSplitCell[i]);
    }
#endif
    
    
    // clean up
    // free memory
    //    printf("a1\n");
    memory->destroy(selCenterLoc);
    //    printf("a2\n");
    memory->destroy(selCenter);
    //    printf("a3\n");
    memory->destroy(selNSceLoc);
    //    printf("a4\n");
    memory->destroy(selNSce);
    //    printf("a5\n");
    memory->destroy(selId);
    //    printf("a6\n");
    memory->destroy(selInd);
    //    printf("a7\n");
    memory->destroy(selR);
    //    printf("a8\n");
    memory->destroy(selRLoc);
    //    printf("a9\n");
    memory->destroy(selEv);
    //    printf("a10\n");
    memory->destroy(selEw);
    //    printf("a11\n");
    memory->destroy(selIndE);
    
}

/* ----------------------------------------------------------------------
 prepare for migration
 - allocate buffer sce's to introduce while removing others
 ---------------------------------------------------------------------- */
void FixSemPM::prepareMig(){

}

/* ----------------------------------------------------------------------
 migration
 - call specified migration procedure
 ---------------------------------------------------------------------- */
void FixSemPM::migrate(){
    migratePolarity();
}

/* ----------------------------------------------------------------------
 migratePolarity
 - select one to take out
 - select one to duplicate
 ---------------------------------------------------------------------- */
void FixSemPM::migratePolarity(){
    // init sth
#ifdef __SEM_DEBUG__
    printf(" - Migrate %d Cells of %d Cells in total\n",nMig,nCell);
#endif
    
    tagMembrane();

    int me;
    MPI_Comm_rank(world,&me);

    int nSel=nMig;
    bool* bSel=bMig;
    
    //++++++++++++++++++
    // GET CELL CENTER
    //++++++++++++++++++
    // allocate stuff
    double **selCenter;
    int * selInd, * selId;
    memory->create(selCenter,nSel,3,"fix/sem:selCenter");
    memory->create(selInd,nSel,"fix/sem:selInd");
    memory->create(selId,nSel,"fix/sem:selId");
    cellCenter(nSel,bSel,selCenter,selId,selInd);	       
    
    //++++++++++++++++++
    // UPDATE POLARITY STATISTICS / THRESHOLDS
    //++++++++++++++++++
    double * selPolThr;
    memory->create(selPolThr,nSel*2,"fix/sem:selPolThr"); 
    updatePolarity(nSel,selId,selInd,selCenter,selPolThr);
    
    //++++++++++++++++++
    // FIND POLY/DEPOLY ELEMENTS
    //++++++++++++++++++
    // allocate stuff
    int * selNSceLoc, * selNSceAll, * selNSceBefore, * selElemLoc;
    memory->create(selNSceLoc,nSel*2,"fix/sem:selNSceLoc");
    memory->create(selNSceAll,nSel*2,"fix/sem:selNSceAll");
    memory->create(selNSceBefore,nSel*2,"fix/sem:selNSceBefore");
    memory->create(selElemLoc,nSel*2,"fix/sem:selElemLoc");

    // init
    for(int s=0;s<nSel;s++) selNSceLoc[s]=selNSceLoc[s+1]=0;
    // count elements above/below threshold, randomly select one
    int nlocal = atom->nlocal;
    double ** x = atom->x;
    int * cell  = atom->cell;
    int *cellpart = atom->cellpart;
    int * mask  = atom->mask;
    for(int i=0;i<nlocal;i++){
	for(int s=0;s<nSel;s++){
	    if(cell[i]==selId[s] && cellpart[i]==0){
		double pol = (selCenter[s][0]-x[i][0])*polarity[selInd[s]][0]+
			     (selCenter[s][1]-x[i][1])*polarity[selInd[s]][1]+
			     (selCenter[s][1]-x[i][2])*polarity[selInd[s]][2];
		if(pol<=selPolThr[s])   selNSceLoc[s]++;
		if(pol>=selPolThr[s+1]) selNSceLoc[s+1]++;
	    }
	}
    }
    MPI_Allreduce(selNSceLoc,selNSceAll,2*nSel,MPI_INT,MPI_SUM,world);
    MPI_Scan(selNSceLoc,selNSceBefore,2*nSel,MPI_INT,MPI_SUM,world);
    // pick random element at root
    if(me==0){
	for(int s=0;s<nSel*2;s++){
	    selElemLoc[s]=rand()%selNSceAll[s]; //!! not yet loc
	}
    }
    //Broadcast it
    MPI_Bcast(selElemLoc, nSel*2, MPI_INT, 0, world);
    
#ifdef __SEM_DEBUG__
    for(int s=0;s<nSel*2;s++){
	printf("%d: SelElemLoc[%d] %d\n",me,s,selElemLoc[s]);
	printf("%d: SelNSceLoc[%d] %d\n",me,s,selNSceLoc[s]);
	printf("%d: SelNSceBefore[%d] %d\n",me,s,selNSceBefore[s]);
	printf("%d: SelNSceAll[%d] %d\n",me,s,selNSceAll[s]);
    }
#endif
    
    
    //++++++++++++++++++
    // find Elements to consider on this CPU
    //++++++++++++++++++
    //allocate

    
    int nSelLoc[2];
    nSelLoc[0]=nSelLoc[1]=0;
    for(int s=0;s<nSel*2;s++){
	if(selElemLoc[s]>=(selNSceBefore[s]-selNSceLoc[s]) && selElemLoc[s]<selNSceBefore[s]){
	    selElemLoc[s]-=(selNSceBefore[s]-selNSceLoc[s]);
	    nSelLoc[s%2]++;
	}else{
	    selElemLoc[s]=-1;
	}
    }
#ifdef __SEM_DEBUG__    
    printf("%d nSelLoc[0] %d, nSelLoc[1] %d\n",me, nSelLoc[0],nSelLoc[1]);
    for(int s=0;s<2*nSel;s++){
	printf("%d selElemLoc[%d] %d\n",me,s,selElemLoc[s]);
    }
#endif
    int cnt0,cnt1;
    cnt0=0;
    cnt1=0;
    
    int * selIndDepolyLoc, * selIndPolyLoc;
    memory->create(selIndDepolyLoc,nSelLoc[0],"fix/sem:selIndDepolyLoc");
    memory->create(selIndPolyLoc,nSelLoc[1],"fix/sem:selIndPolyLoc");
    for(int s=0;s<nSel;s++){
	if(selElemLoc[s]>=0){
	    selIndDepolyLoc[cnt0]=s;
	    cnt0++;
	}
	if(selElemLoc[s+1]>=0){
	    selIndPolyLoc[cnt1]=s;
	    cnt1++;
	}
    }

#ifdef __SEM_DEBUG__
    for(int s=0;s<nSelLoc[0];s++){
	printf("%d, Deoly selIndDeployLoc[%d] %d \n",me,s,selIndDepolyLoc[s]);
    }
    for(int s=0;s<nSelLoc[1];s++){
	printf("%d, Poly selIndPolyLoc[%d] %d\n",me,s,selIndPolyLoc[s]);
    }
#endif
    
    //++++++++++++++++++
    // HANDLE POLY
    // create new particles locally
    //++++++++++++++++++

    bigint natoms_previous = atom->natoms;
    int nlocal_previous = atom->nlocal;
    
    if(nSelLoc[1]>0){
	int cnt=0;
	for(int i=0;i<nlocal;i++){
	    for(int l=0;l<nSelLoc[1];l++){
		int s=selIndPolyLoc[l];
		int c=selInd[s];
		//printf("%d: create: l: %d,s: %d,c %d\n",me,l,s,c);
		if(cell[i]==selId[s] && cellpart[i]==0){
		    double pol = (selCenter[s][0]-x[i][0])*polarity[c][0]+
				 (selCenter[s][1]-x[i][1])*polarity[c][1]+
				 (selCenter[s][2]-x[i][2])*polarity[c][2];
		    if(pol>=selPolThr[2*s+1]){
			if(selElemLoc[2*s+1]==0){
//			    printf("CREATED\n");
			    //++++++++++++++
			    // create element
			    //++++++++++++++
			    double coord[3];
			    coord[0]=x[i][0];
			    coord[1]=x[i][1];
			    coord[2]=x[i][2];
			
			    atom->avec->create_atom(tyCell[c],coord);
			    int m=atom->nlocal-1;
			    atom->cell[m]=idCell[c];
			    atom->p[m]=0;
			    atom->mask[m]|=polyMGroupBitmask;
			    //++++++++++++++
			    // done create element
			    //++++++++++++++
			    selElemLoc[2*s+1]--;
			    cnt++;
			    if(cnt==nSelLoc[1]) break;
			}else{
			    selElemLoc[2*s+1]--;
			}
		    }
		}
	    }
	    if(cnt==nSelLoc[1]) break;
	}
	// call some fixes
	for (int m = 0; m < modify->nfix; m++) {
	    Fix *fix = modify->fix[m];
	    if (fix->create_attribute){
		for (int i = nlocal_previous; i < nlocal; i++){
		    printf("EVER CALLED?");
		    fix->set_arrays(i);
		}
	    }
	}
    }    
            
    //++++++++++++++++++
    // HANDLE DEPOLY
    // tag them
    //++++++++++++++++++

    if(nSelLoc[0]>0){
	int cnt=0;
	for(int i=0;i<nlocal;i++){
	    for(int l=0;l<nSelLoc[0];l++){
		int s=selIndDepolyLoc[l];
		int c=selInd[s];
		//printf("%d: takeout: l: %d,s: %d,c %d\n",me,l,s,c);
		if(cell[i]==selId[s] && cellpart[i]==0){
		    double pol = (selCenter[s][0]-x[i][0])*polarity[c][0]+
				 (selCenter[s][1]-x[i][1])*polarity[c][1]+
				 (selCenter[s][2]-x[i][2])*polarity[c][2];
		    if(pol<=selPolThr[2*s]){
			if(selElemLoc[2*s]==0){
//			    printf("TAGGED OUT\n");
			    //++++++++++++++
			    // tag element to go out
			    // remove other masks (except for all)
			    //++++++++++++++
			    mask[i]=1;
			    mask[i]|=depolyMGroupBitmask;
			    //++++++++++++++
			    // done tag
			    //++++++++++++++
			    selElemLoc[2*s]--;
			    cnt++;
			    if(cnt==nSelLoc[0]) break;
			}else{
			    selElemLoc[2*s]--;
			}
		    }
		}
	    }
	    if(cnt==nSelLoc[0]) break;
	}
    }
    // new total # of atoms
    bigint rlocal = atom->nlocal;
    MPI_Allreduce(&rlocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
    if (atom->natoms < 0 || atom->natoms > MAXBIGINT)
	error->all(FLERR,"Too many total atoms");
    
//    printf("atoms: %d\n",atom->natoms);
    
    // set tag # of new particles beyond all previous atoms
    // if global map exists, reset it now instead of waiting for comm
    // since deleting atoms messes up ghosts
    if (atom->tag_enable) {
      atom->tag_extend();
      if (atom->map_style) {
	atom->nghost = 0;
	atom->map_init();
	atom->map_set();
      }
    }
    
    // never invoced
    if (atom->molecular) {
	int **nspecial = atom->nspecial;
	for (int i = nlocal_previous; i < atom->nlocal; i++) {
	    nspecial[i][0] = 0;
	    nspecial[i][1] = 0;
	    nspecial[i][2] = 0;
	}
    }
    
    //@@@SEM: force rebuild of neighbour List next time step:
    //needs neighbour must_check enabled!!!
    nPolymerizeSce+=nSel*2;
    next_reneighbor=update->ntimestep+1;
    
    
    //++++++++++++++++++
    // DEALLOCATE THINGS
    //++++++++++++++++++
    memory->destroy(selCenter);
    memory->destroy(selInd);
    memory->destroy(selId);
    memory->destroy(selPolThr);
    memory->destroy(selNSceLoc);
    memory->destroy(selNSceAll);
    memory->destroy(selNSceBefore);
    memory->destroy(selElemLoc);
    memory->destroy(selIndPolyLoc);
    memory->destroy(selIndDepolyLoc);    
}
/* ----------------------------------------------------------------------
 tag membrane
 ---------------------------------------------------------------------- */
void FixSemPM::tagMembrane(){
    //tagMembraneDens();
}

/* ----------------------------------------------------------------------
 update polarity
 ---------------------------------------------------------------------- */
void FixSemPM::updatePolarity(){
    // UPDATE POLARITY VECTOR
    // - right now leave constant
    
    // UPDATE POLARITY OF ELEMENTS
    // - determine cell center
    // - determine polarity as projection of distance elem-cell center
    //   in direction of polarity vector
}

/* ----------------------------------------------------------------------
    tag membrane desnity based
    - run through elements, if Dens<11, mark as membrane
---------------------------------------------------------------------- */
void FixSemPM::tagMembraneDens(){
    int nlocal = atom->nlocal;
    double * dens = atom->dens;
    int * cellpart = atom->cellpart;
    for(int i=1;i<nlocal;i++){
	cellpart[i]=(dens[i]<11)?1:0;
    }
}


/* ----------------------------------------------------------------------
-------------------------------------------------------------------------
    HELPER FUNCTIONS
-------------------------------------------------------------------------    
 ---------------------------------------------------------------------- */
void FixSemPM::updatePolarity(int nSel, int *selId, int* selInd, double** selCenter,double *selPolThr){
    
    // stuff
    int me;
    MPI_Comm_rank(world,&me);
    
    // allocate
    double * selPolStat, * selPolStatLoc, * selPolMin, * selPolMinLoc, * selPolMax, * selPolMaxLoc;
    memory->create(selPolStat,nSel*2,"fix/sem:selPolStat");
    memory->create(selPolStatLoc,nSel*2,"fix/sem:selPolStatLoc");
    memory->create(selPolMin,nSel,"fix/sem:selPolMin");
    memory->create(selPolMinLoc,nSel,"fix/sem:selPolMinLoc");
    memory->create(selPolMax,nSel,"fix/sem:selPolMax");
    memory->create(selPolMaxLoc,nSel,"fix/sem:selPolMaxLoc");
    
    // init
    for(int s=0;s<nSel;s++){
	selPolMinLoc[s]=0.0;
	selPolMaxLoc[s]=0.0;
	selPolStatLoc[s]=0.0;
	selPolStatLoc[s+1]=0.0;
    }
    //accumulate locally
    int nlocal = atom->nlocal;
    int *cell  = atom->cell;
    int *cellpart = atom->cellpart;
    double **x = atom->x;    
    for(int i=0;i<nlocal;i++){
	for(int s=0;s<nSel;s++){
	    if(cell[i]==selId[s] && cellpart[i]==0){
		double pol =(selCenter[s][0]-x[i][0])*polarity[selInd[s]][0]+
			    (selCenter[s][1]-x[i][1])*polarity[selInd[s]][1]+
			    (selCenter[s][2]-x[i][2])*polarity[selInd[s]][2];
		selPolStatLoc[s]+=pol;
		selPolStatLoc[s+1]+=pol*pol;
		selPolMinLoc[s]=(selPolMinLoc[s]<pol)?selPolMinLoc[s]:pol;
		selPolMaxLoc[s]=(selPolMaxLoc[s]>pol)?selPolMaxLoc[s]:pol;
	    }
	}
    }
    // accumulate statistics globally
    MPI_Allreduce(selPolStatLoc,selPolStat,nSel*2,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(selPolMinLoc,selPolMin,nSel,MPI_DOUBLE,MPI_MIN,world);
    MPI_Allreduce(selPolMaxLoc,selPolMax,nSel,MPI_DOUBLE,MPI_MAX,world);
    for(int s=0;s<nSel;s++){
	selPolStat[s]/=nSceCell[selInd[s]];
	selPolStat[s+1]/=nSceCell[selInd[s]];
	selPolStat[s+1]-=selPolStat[s]*selPolStat[s];
    }
    
    //++++++++++++++++++
    // FIND THRESHOLD LEVELS FOR POLY/DEPOLY
    //++++++++++++++++++
    // choose threshold level for poly/depoly based on statistics
    // do at root
    double DepolyFac=0.2;
    double PolyFac=0.2;   
    if(me==0){
	for(int s=0;s<nSel;s++){
	    selPolThr[s]=MAX(selPolMin[s],selPolStat[s]-selPolStat[s+1]*DepolyFac); 
	    selPolThr[s+1]=MIN(selPolMax[s],selPolStat[s]+selPolStat[s+1]*PolyFac);
	}
    }
    // distribute thresholds;
    MPI_Bcast(selPolThr,nSel*2,MPI_DOUBLE,0,world);
#ifdef __SEM_DEBUG__
    for(int s=0;s<2*nSel;s++){
	printf("Polarity[%d] %f\n",s,selPolThr[s]);
    }
#endif
    // free memory
    memory->destroy(selPolStat);
    memory->destroy(selPolStatLoc);
    memory->destroy(selPolMin);
    memory->destroy(selPolMinLoc);
    memory->destroy(selPolMax);
    memory->destroy(selPolMaxLoc);
}

void FixSemPM::cellCenter(int nSel, bool *bSel, double ** selCenter, int *selId,int *selInd){

    double **selCenterLoc;
    memory->create(selCenterLoc,nSel,3,"fix/sem:selCenterLoc");
    
    int cur=0;
    for(int c=0;c<nCell;c++){
	if(bSel[c]){
	    selCenterLoc[cur][0]=0.0;
	    selCenterLoc[cur][1]=0.0;
	    selCenterLoc[cur][2]=0.0;
	    selId[cur]=idCell[c];
	    selInd[cur]=c;
	    cur++;
	}
    }
    
    // find centers of migrating 
    // accumulate locations localy
    int nlocal = atom->nlocal;
    int *cell = atom->cell;
    int *mask  = atom->mask;
    double **x = atom->x;
    
    for (int i=0;i<nlocal;i++){	
	for(int j=0;j<nSel;j++){
	    if(cell[i]==selId[j]){ // @@@??? maybe just do for all
		selCenterLoc[j][0]+=x[i][0];
		selCenterLoc[j][1]+=x[i][1];
		selCenterLoc[j][2]+=x[i][2];
	    }
	}
    }
    
    MPI_Allreduce(&selCenterLoc[0][0],&selCenter[0][0],3*nSel,MPI_DOUBLE,MPI_SUM,world);
    
    // scale
    for(int i=0;i<nSel;i++){
	//	printf("proc: %d cell: %d | center: %f %f %f\n",me,selId[i],selCenter[i][0],selCenter[i][1],selCenter[i][2]);
	selCenter[i][0]/=(double)nSceCell[selInd[i]];
	selCenter[i][1]/=(double)nSceCell[selInd[i]];
	selCenter[i][2]/=(double)nSceCell[selInd[i]];
    }
    memory->destroy(selCenterLoc);
    
}

