/*
 *  pair_semNucM.cpp
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 1/26/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include "pair_semNucM.h"

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

#include "update.h"
using namespace LAMMPS_NS;

#define EPS 1e-5

/* ---------------------------------------------------------------------- */

PairSemNucM::PairSemNucM(LAMMPS *lmp) : Pair(lmp) {}

/* ----------------------------------------------------------------------
 // what was allocated
 
 // what to allocate
 
 // genral
 setflag[ii,jj]
 cutsq[ii,jj]
 
 r0[cellType][NucCellType];
 cutDens2[cellType][NucCellType];
 
 u0[cellType][NucCellType]
 rho[cellType][NucCellType]
 alpha[cellType][NucCellType]
 rShift[cellType][NucCellType]
 rhoScl[cellType][NucCellType]
 semp1[cellType][NucCellTyp
 rShiftSq[cellType][NucCellType]
 rMinSq[cellType][NucCellType]
 
 
 --------------------------------------------------------------------------*/

PairSemNucM::~PairSemNucM()
{
    if (allocated) {
	// general
        memory->destroy(setflag);
        memory->destroy(cutsq);
        
        memory->destroy(r0);
	
	// dens
	memory->destroy(cutDens2);
	
	// intra
	memory->destroy(u0);
	memory->destroy(rho);
	memory->destroy(alpha);
	memory->destroy(rShift);	
	memory->destroy(rhoScl);
	memory->destroy(semp1);
	memory->destroy(rShiftSq);
	memory->destroy(rMinSq);
	memory->destroy(cut);
	
	// NO INTER
	
	// ??
	memory->destroy(offset);
    }
}

/* ---------------------------------------------------------------------- */

void PairSemNucM::compute(int eflag, int vflag)
{
    int i,j,ii,jj,inum,jnum,itype,jtype,icell,jcell;
    double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair,ip,jp;
    double rsq,rsqSclP,pFac,E,factor_lj;
    int *ilist,*jlist,*numneigh,**firstneigh;
    
    evdwl = 0.0;
    if (eflag || vflag) ev_setup(eflag,vflag);
    else evflag = vflag_fdotr = 0;
    
    double **x = atom->x;
    double **f = atom->f;
    int *type = atom->type;
    int *cell  = atom->cell;
    double *p  = atom->p;
    double *dens = atom->dens;
    int nlocal = atom->nlocal;
    int nall = nlocal + atom->nghost;
    double *special_lj = force->special_lj;
    int newton_pair = force->newton_pair;
    
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    
    // loop over neighbors of my atoms
#define __SEM_DEBUG__NOT
#ifdef __SEM_DEBUG__
    int me;
    MPI_Comm_rank(world,&me);
#endif
    
    for (ii = 0; ii < inum; ii++) {
	
	i = ilist[ii];
	xtmp = x[i][0];
	ytmp = x[i][1];
	ztmp = x[i][2];
	itype = type[i];
	icell  = cell[i];
	ip     = p[i];
	//dens[i]=0.0;	//@@@??? do somewhere else if newton is on - for now in fix_proliferation
#ifdef __SEM_DEBUG__
	printf("++i" BIGINT_FORMAT " | proc: %d, particle %d, location: %f %f %f\n",update->ntimestep, me, atom->tag[i],x[i][0],x[i][1],x[i][2]);
#endif
	
	jlist = firstneigh[i];
	jnum = numneigh[i];
	//	printf("touch i: %i - in list: %i \n",i,jnum);
	for (jj = 0; jj < jnum; jj++) {
#ifdef __SEM_DEBUG__
	    printf("++j" BIGINT_FORMAT " | proc: %d, particle %d, location: %f %f %f\n",update->ntimestep, me, atom->tag[i],x[i][0],x[i][1],x[i][2]);
#endif
	    j = jlist[jj];	
	    if (j < nall) factor_lj = 1.0;
	    else {
		factor_lj = special_lj[j/nall];
		j %= nall;
	    }
	    
	    delx = xtmp - x[j][0];
	    dely = ytmp - x[j][1];
	    delz = ztmp - x[j][2];
	    jtype = type[j];
	    jcell = cell[j];
	    jp    = p[j];
	    rsq = (delx*delx+dely*dely+delz*delz);
	    if (rsq < cutsq[itype][jtype]) {
		// intra
		if(icell == jcell){
		    // only consider p scaling for nuc-nuc .. 
		    //pFac=(itype==jtype)?MIN(ip,jp)+EPS:1.0;
		    pFac=MIN(ip,jp)+EPS;
		    double rs2=sqrt(rsq);
		    rs2-=rShift[itype][jtype]*pFac;
		    rs2*=rs2;

		    rsqSclP=MAX(rMinSq[itype][jtype],rs2)/(pFac*pFac);    

		    E = exp(rho[itype][jtype]-rhoScl[itype][jtype]*rsqSclP);
		    fpair = factor_lj * semp1[itype][jtype]*(2.*E*E-alpha[itype][jtype]*E) * pFac;
		    f[i][0] += delx*fpair; // delx, dely, delz is scaled with last factore   <                >
		    f[i][1] += dely*fpair;
		    f[i][2] += delz*fpair;
		    
		    // accumulate dens
		    if(rsq<cutDens2[itype][jtype])dens[i]+=3.0;
		    
		    if (newton_pair || j < nlocal) {
			f[j][0] -= delx*fpair;
			f[j][1] -= dely*fpair;
			f[j][2] -= delz*fpair;
			//if(rsq<cutDens2[itype])
			if(rsq<cutDens2[itype][jtype])dens[j]+=3.0;
		    }
		    
		    if (eflag) {
			evdwl = u0[itype][jtype] * (E*E - alpha[itype][jtype]*E) -
			offset[itype][jtype];
			evdwl *= factor_lj;
		    }
		    
		}
		if (evflag) ev_tally(i,j,nlocal,newton_pair,
				     evdwl,0.0,fpair,delx,dely,delz);
	    }
	}
    }
    
    if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
 allocate all arrays 
 
 // what to allocate
 
 // genral
 setflag[ii,jj]
 cutsq[ii,jj]
 
 r0[cellType][NucCellType];
 cutDens2[cellType][NucCellType];
 
 // intra
 u0[cellType][NucCellType]
 rho[cellType][NucCellType]
 alpha[cellType][NucCellType]
 rShift[cellType][NucCellType]
 rhoScl[cellType][NucCellType]
 semp1[cellType][NucCellTyp
 rShiftSq[cellType][NucCellType]
 rMinSq[cellType][NucCellType]
 
 ------------------------------------------------------------------------- */

void PairSemNucM::allocate()
{
    allocated = 1;
    int n = atom->ntypes; // cell types
    
    // general
    memory->create(setflag,n+1,n+1,"pair:setflag");
    for (int i = 1; i <= n; i++)
	for (int j = i; j <= n; j++)
	    setflag[i][j] = 0;
    
    memory->create(cutsq,n+1,n+1,"pair:cutsq");
    memory->create(r0,n+1,n+1,"pair:r0");
    
    // dens
    memory->create(cutDens2,n+1,n+1,"pair::cutDens2");
    
    //intra
    memory->create(rho,n+1,n+1,"pair:rho");
    memory->create(u0,n+1,n+1,"pair:u0");
    memory->create(alpha,n+1,n+1,"pair:alpha");
    memory->create(rShift,n+1,n+1,"pair:rShift");
    memory->create(rhoScl,n+1,n+1,"pair:rhoScl");
    memory->create(semp1,n+1,n+1,"pair:semp1");
    memory->create(rShiftSq,n+1,n+1,"pair:rShiftSq");
    memory->create(rMinSq,n+1,n+1,"pair:rMinSq");
    memory->create(cut,n+1,n+1,"pair:cut");
        
    //??
    memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
 global settings 
 ------------------------------------------------------------------------- */

void PairSemNucM::settings(int narg, char **arg)
{
    if (narg != 1) error->all(FLERR,"Illegal pair_style command");
    
    cut_global = utils::numeric(FLERR,arg[0],false,lmp);
    
    // reset cutoffs that have been explicitly set
    
    if (allocated) {
	for (int i = 1; i <= atom->ntypes; i++){
	    for (int j = i+1; j <= atom->ntypes; j++)
		if (setflag[i][j]) cut[i][j] = cut_global;
	    
	}
    }
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 nargs: 6(7)
 args: [inter/intra to differentiate between inter/intra coefficients]
 -0) atom_type	(CELL TYPE)
 -1) atom_type  (CELL TYPE NUCLEUS)
 -2) u0 	(u_0)
 -3) rho	[internally, store alpha=rho/req2]
 -4) alpha
 -5) r0			(r_eq)[store r0^2]
 -6) densFac		(fraction of r0)
 -7) rShift             (shift of potential, store shift^2)
 -8) cutoff		(optional)
 we store internally
 rho:	rho
 u0:		u0
 alpha:  alpha
 rhoScl: rho/r0^2
 semp1:  2*u0*rho/r0^2
 cut:	cutoff
 r0:				r_eq
 
 ------------------------------------------------------------------------- */

void PairSemNucM::coeff(int narg, char **arg)
{
    int ii,jj;
    if (narg < 8 || narg > 9) error->all(FLERR,"Incorrect args for pair coefficients a");
    if (!allocated) allocate();
    
    int ilo,ihi,jlo,jhi;
    utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
    utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);
        
    double u0_one = utils::numeric(FLERR,arg[2],false,lmp);
    double rho_one = utils::numeric(FLERR,arg[3],false,lmp);
    double alpha_one = utils::numeric(FLERR,arg[4],false,lmp);
    double r0_one = utils::numeric(FLERR,arg[5],false,lmp);
    double densFac_one = utils::numeric(FLERR,arg[6],false,lmp);
    double rShift_one = utils::numeric(FLERR,arg[7],false,lmp);
    
    double cut_one = cut_global;
    if (narg == 9) cut_one = utils::numeric(FLERR,arg[8],false,lmp);
    
    
    int count = 0;
    for (int i=ilo; i<=ihi; i++){
	for (int j= MAX(jlo,i); j<=jhi; j++){
	    // set general coefficients
	    setflag[i][j]=1;
	    r0[i][j]=r0_one;
	    // dens
	    cutDens2[i][j]=r0_one*r0_one*densFac_one*densFac_one;
    
	    // set intra coefficients	
	    u0[i][j]=u0_one;
	    rho[i][j]=rho_one;
	    alpha[i][j]=alpha_one;
	    rShift[i][j]=rShift_one;
	    cut[i][j]=cut_one;
	    count++;
	    printf("** pair_semNucM: set coeffs: u0 %f, rho0 %f, alpha %f, r0 %f, shift %f, cut %f\n",
		   u0[i][j],rho[i][j],alpha[i][j],r0[i][j],rShift[i][j],cut[i][j]);
	    
	}
    }
    if (count==0) error->all(FLERR,"Incorrect args for pair coefficients b");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSemNucM::init_one(int i, int j)
{	   
    printf("** pair_semNucM: init_one for %d and %d\n",i,j);
    double r02;
    if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
    r02=r0[i][j]*r0[i][j];
    rhoScl[i][j]=rho[i][j]/r02;
    semp1[i][j]=2.*rho[i][j]*u0[i][j]/r02;
    rShiftSq[i][j]=rShift[i][j]*rShift[i][j];
    rMinSq[i][j]=(rShift[i][j]>0)?0.01*r02:0.0;
    
    offset[i][j]=0.0; //@@@
    
    //copy stuff
    cutDens2[j][i]=cutDens2[i][j];
    u0[j][i]=u0[i][j];
    rho[j][i]=rho[i][j];
    alpha[j][i]=alpha[i][j];
    rShift[j][i]=rShift[i][j];
    cut[j][i]=cut[i][j];
    rhoScl[j][i]=rhoScl[i][j];
    semp1[j][i]=semp1[i][j];
    rShiftSq[j][i]=rShiftSq[i][j];
    rMinSq[j][i]=rMinSq[i][j];
    offset[j][i]=offset[i][j];

    return cut[i][j];
}

/* ----------------------------------------------------------------------
 proc 0 writes to restart file
 ------------------------------------------------------------------------- */
//@@@SEM RESTART NOT TESTED
void PairSemNucM::write_restart(FILE *fp)
{
    write_restart_settings(fp);
    
    int i,j;
    for(i=1; i<=atom->ntypes; i++){
	for(j=1; j<atom->ntypes; j++){
	    fwrite(&setflag[i][j],sizeof(int),1,fp);
	    if(setflag[i][j]){
		//intra
		fwrite(&u0[i][j],sizeof(double),1,fp);
		fwrite(&rho[i][j],sizeof(double),1,fp);
		fwrite(&alpha[i][j],sizeof(double),1,fp);
		fwrite(&rShift[i][j],sizeof(double),1,fp);
		fwrite(&cut[i][j],sizeof(double),1,fp);	    
	    }
	}
    }
}

/* ----------------------------------------------------------------------
 proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */
//@@@SEM RESTART NOT TESTED
void PairSemNucM::read_restart(FILE *fp)
{
    read_restart_settings(fp);
    
    allocate();
    
    int i,j;
    int me = comm->me;
    for (i=1; i<=atom->ntypes; i++){
	for (j=i; j<=atom->ntypes; j++){
	    if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
	    MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
	    if (setflag[i][j]) {
		if (me == 0){
		    fread(&u0[i][j],sizeof(double),1,fp);
		    fread(&rho[i][j],sizeof(double),1,fp);
		    fread(&alpha[i][j],sizeof(double),1,fp);
		    fread(&rShift[i][j],sizeof(double),1,fp);
		    fread(&cut[i][j],sizeof(double),1,fp);
		}
		MPI_Bcast(&u0[i][j],1,MPI_DOUBLE,0,world);
		MPI_Bcast(&rho[i][j],1,MPI_DOUBLE,0,world);
		MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
		MPI_Bcast(&rShift[i][j],1,MPI_DOUBLE,0,world);
		MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
	    }
	}
    }
}

/* ----------------------------------------------------------------------
 proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void PairSemNucM::write_restart_settings(FILE *fp)
{	
    fwrite(&cut_global,sizeof(double),1,fp);
    fwrite(&offset_flag,sizeof(int),1,fp);
    fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
 proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void PairSemNucM::read_restart_settings(FILE *fp)
{
    if (comm->me == 0) {
	fread(&cut_global,sizeof(double),1,fp);
	fread(&offset_flag,sizeof(int),1,fp);
	fread(&mix_flag,sizeof(int),1,fp);
    }
    MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
    MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */
//@@@SEM where is this called?
//how to pass the polymerization factor and the cell information?
double PairSemNucM::single(int i, int j, int itype, int jtype, int icell, int jcell,double poly, double rsq,
			  double factor_coul, double factor_lj,
			  double &fforce)
{
    double r,E,phi;
    
    ////////
    printf("ATTENTION _ DONE!/n");
    rsq=MAX(rMinSq[itype][jtype],rsq-rShiftSq[itype][jtype]);
    r=sqrt(rsq);
    rsq/=(poly*poly+EPS);
    
    if(icell == jcell){
	E = exp(rho[itype][jtype]-rhoScl[itype][jtype]*rsq);
	fforce = factor_lj * semp1[itype][jtype]*(2.*E*E - alpha[itype][jtype]*E)/r;
	phi=u0[itype][jtype] * (E*E - alpha[itype][jtype]*E) - offset[itype][jtype];	
    }
    return factor_lj*phi;
}
