/*
 *  pari_semp.cpp
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 12/9/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_semp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSemp::PairSemp(LAMMPS *lmp) : Pair(lmp) {}

/* ----------------------------------------------------------------------
// what was allocated

 // what to allocate
 
 // genral
 setflag[ii,jj]
 cut[ii,jj]
 cutsq[ii,jj]
 r0[ii,jj];
 
 // intra
 cut_loc[cellType]
 r02[cellType]
 rhoScl[cellType]
 u0[cellType]
 semp1[cellType]
 
 
 // inter tmp
 u0_inter[ii,jj]
 rho_inter[cellType]
 r02_inter[ii,jj]
 rhoScl_inter[ii,jj]
 semp1_inter[ii,jj

//??
morse1
offset
--------------------------------------------------------------------------*/

PairSemp::~PairSemp()
{
	if (allocated) {
		std::cout<<"destructor"<<std::endl;
		// general
		memory->destroy(setflag);
		memory->destroy(cut);
		memory->destroy(cutsq);
		memory->destroy(r0);
		
		// intra
		memory->destroy(cut_loc);
		memory->destroy(rhoScl);
		memory->destroy(u0);
		memory->destroy(semp1);
		
		// inter
		memory->destroy(u0_inter);
		memory->destroy(rho_inter);
		memory->destroy(r02_inter);
		memory->destroy(rhoScl_inter);
		memory->destroy(semp1_inter);
		
		// ??
		memory->destroy(offset);
	}
}

/* ---------------------------------------------------------------------- */

void PairSemp::compute(int eflag, int vflag)
{
//	std::cout<<"compute"<<std::endl;
	int i,j,ii,jj,inum,jnum,itype,jtype,imol,jmol;
	double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
	double rsq,r,dr,dexp,factor_lj;
	int *ilist,*jlist,*numneigh,**firstneigh;
	
	evdwl = 0.0;
	if (eflag || vflag) ev_setup(eflag,vflag);
	else evflag = vflag_fdotr = 0;
	
	double **x = atom->x;
	double **f = atom->f;
	int *type = atom->type;
	int *mol  = atom->molecule;//@@@
	int nlocal = atom->nlocal;
	int nall = nlocal + atom->nghost;
	double *special_lj = force->special_lj;
	int newton_pair = force->newton_pair;
	
	inum = list->inum;
	ilist = list->ilist;
	numneigh = list->numneigh;
	firstneigh = list->firstneigh;
	
	// loop over neighbors of my atoms
	
	for (ii = 0; ii < inum; ii++) {
		i = ilist[ii];
		xtmp = x[i][0];
		ytmp = x[i][1];
		ztmp = x[i][2];
		itype = type[i];
		imol  = mol[i]; //@@@
		jlist = firstneigh[i];
		jnum = numneigh[i];
		
		for (jj = 0; jj < jnum; jj++) {
			j = jlist[jj];
			
			if (j < nall) factor_lj = 1.0;
			else {
				factor_lj = special_lj[j/nall];
				j %= nall;
			}
			
			delx = xtmp - x[j][0];
			dely = ytmp - x[j][1];
			delz = ztmp - x[j][2];
			rsq = delx*delx + dely*dely + delz*delz;
			jtype = type[j];
			jmol = mol[j];//@@@
			
			if (rsq < cutsq[itype][jtype]) {
				dr = rsq - r02_inter[itype][jtype];
				//r=sqrt(rsq);//@@@
				//dr = r - r0[itype][jtype];
				// inter
				if(imol == jmol){
					dexp = exp(-rhoScl[itype] * dr);
					fpair = factor_lj * semp1[itype]*(dexp*dexp - dexp);
					//fpair = factor_lj * semp1[itype]*(dexp*dexp - dexp)/(r*2);

					f[i][0] += delx*fpair;
					f[i][1] += dely*fpair;
					f[i][2] += delz*fpair;
					if (newton_pair || j < nlocal) {
						f[j][0] -= delx*fpair;
						f[j][1] -= dely*fpair;
						f[j][2] -= delz*fpair;
					}
					
					if (eflag) {
						evdwl = u0[itype] * (dexp*dexp - 2.0*dexp) -
						offset[itype][jtype];
						evdwl *= factor_lj;
					}
					
					
				}else{
					dexp = exp(-rhoScl_inter[itype][jtype]*dr);
					fpair = factor_lj * semp1_inter[itype][jtype]*(dexp*dexp - dexp);
					//fpair = factor_lj * semp1_inter[itype][jtype]*(dexp*dexp - dexp)/(r*2);
					f[i][0] += delx*fpair;
					f[i][1] += dely*fpair;
					f[i][2] += delz*fpair;
					if (newton_pair || j < nlocal) {
						f[j][0] -= delx*fpair;
						f[j][1] -= dely*fpair;
						f[j][2] -= delz*fpair;
					}
					if (eflag) {
						evdwl = u0_inter[itype][itype] * (dexp*dexp - 2.0*dexp) -
						offset[itype][itype];
						evdwl *= factor_lj;
					}
										
					
				}
				if (evflag) ev_tally(i,j,nlocal,newton_pair,
									 evdwl,0.0,fpair,delx,dely,delz);
					
				
			}
		}
	}
	
	if (vflag_fdotr) virial_fdotr_compute();
//	std::cout<<"DONE compute"<<std::endl;
}

/* ----------------------------------------------------------------------
 allocate all arrays 
 
 // what to allocate
 
 // genral
 setflag[ii,jj]
 cut[ii,jj]
 cutsq[ii,jj]
 r0[ii,jj]
 
 // intra
 cut_loc[cellType]
 r02[cellType]
 rhoScl[cellType]
 u0[cellType]
 semp1[cellType]

 
 // inter tmp
 u0_inter[ii,jj]
 rho_inter[cellType]
 r02_inter[ii,jj]
 rhoScl_inter[ii,jj]
 semp1_inter[ii,jj]
 
 //??
 offset
 ------------------------------------------------------------------------- */

void PairSemp::allocate()
{
	std::cout<<"allocate"<<std::endl;
	allocated = 1;
	int n = atom->ntypes; // cell types
	
	// general
        memory->create(setflag,n+1,n+1,"pair:setflag");
	for (int i = 1; i <= n; i++)
		for (int j = i; j <= n; j++)
			setflag[i][j] = 0;
	
	memory->create(cut,n+1,n+1,"pair:cut");
        memory->create(cutsq,n+1,n+1,"pair:cutsq");
        memory->create(r0,n+1,n+1,"pair:r0");
	
	//intra
	memory->create(cut_loc,n,"pair:cut_loc");
	memory->create(r02,n,"pair:r02");
        memory->create(rhoScl,n,"pair:rhoScl");
        memory->create(u0,n,"pair:u0");
	memory->create(semp1,n,"pair:semp1");
	
	//intra tmp
        memory->create(u0_inter,n+1,n+1,"pair:u0_inter");
        memory->create(rho_inter,n,"pair:rho_inter");
        memory->create(r02_inter,n+1,n+1,"pair:r02_inter");
        memory->create(rhoScl_inter,n+1,n+1,"pair:rhoScl_inter");
        memory->create(semp1_inter,n+1,n+1,"pair:semp1_inter");
	
	//??
        memory->create(offset,n+1,n+1,"pair:offset");
	std::cout<<"DONE allocate"<<std::endl;
}

/* ----------------------------------------------------------------------
 global settings 
 ------------------------------------------------------------------------- */

void PairSemp::settings(int narg, char **arg)
{
	std::cout<<"settings"<<std::endl;
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
	-1) ignored	(just to be coherent with other atom styles)
	-1) u0_intra	(u_0)
	-2) u0_inter
	-3) rho_intra	[internally, store alpha=rho/req2]
	-4) rho_inter
	-5) r0			(r_eq)[store r0^2]
	-6) cutoff		(optional)
we store intenally
	U0_intra=4*rho*u0/req2
	RHO_intra=rho/req2
	r02=req2

 ------------------------------------------------------------------------- */

void PairSemp::coeff(int narg, char **arg)
{
	std::cout<<"coeff"<<std::endl;	
	int ii,jj;
	if (narg < 7 || narg > 8) error->all(FLERR,"Incorrect args for pair coefficients");
	if (!allocated) allocate();
	
	int cellType=utils::inumeric(FLERR,arg[0],false,lmp);
	int spare = utils::inumeric(FLERR,arg[1],false,lmp);
	printf("PAIRCOEFF_SEM: second parameter %d ignored.. mixing rules applied\n",spare);
	double u0_intra_one = utils::numeric(FLERR,arg[2],false,lmp);
	double u0_inter_one = utils::numeric(FLERR,arg[3],false,lmp);
	double rho_intra_one = utils::numeric(FLERR,arg[4],false,lmp);
	double rho_inter_one = utils::numeric(FLERR,arg[5],false,lmp);
	double r0_one = utils::numeric(FLERR,arg[6],false,lmp);
	
	double cut_one = cut_global;
	if (narg == 8) cut_one = utils::numeric(FLERR,arg[7],false,lmp);

	// set general coefficients
	setflag[cellType][cellType]=1;
	r0[cellType][cellType]=r0_one;
	
	// set intra coefficients
	cut_loc[cellType]=cut_one;
	r02[cellType]=r0_one*r0_one; //@@@ obsolete?
	rhoScl[cellType]=rho_intra_one/r02[cellType];
	u0[cellType]=u0_intra_one;
	semp1[cellType]=4.*rho_intra_one*u0_intra_one/r02[cellType];
	
	// store coefficients to interpolate inter coefficients
	u0_inter[cellType][cellType]=u0_inter_one;
	rho_inter[cellType]=rho_inter_one;

	// update intra coefficient table
	for (int i=1;i<=atom->ntypes;i++) {
		// only if other potential is set
		if(setflag[i][i]) {
			setflag[cellType][i]=1;
			ii=MIN(i,cellType);
			jj=MAX(i,cellType);
			u0_inter[ii][jj]=0.5*(u0_inter[ii][jj]+u0_inter[ii][jj]);
			r0[ii][jj]=0.5*(r0[ii][ii]+r0[jj][jj]);
			r02_inter[ii][jj]=0.5*(r0[ii][ii]*r0[ii][ii]+r0[jj][jj]*r0[jj][jj]);
			rhoScl_inter[ii][jj]=0.5*(rho_inter[ii]+rho_inter[jj])/(r02_inter[ii][jj]);
			semp1_inter[ii][jj]=4.0*u0_inter[ii][jj]*rhoScl_inter[ii][jj];
			cut[ii][jj]=MAX(cut_loc[ii],cut_loc[jj]);
		}
	};
	std::cout<<"DONE coeff"<<std::endl;
}


/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSemp::init_one(int i, int j)
{
	std::cout<<"init_one"<<std::endl;
	
	if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
	
	//morse1[i][j] = 2.0*d0[i][j]*alpha[i][j];
	
//	if (offset_flag) {
//		double alpha_dr = -alpha[i][j] * (cut[i][j] - r0[i][j]);
//		offset[i][j] = d0[i][j] * (exp(2.0*alpha_dr) - 2.0*exp(alpha_dr));
//	} else offset[i][j] = 0.0;
	offset[i][j] = 0.0;
	
	u0_inter[j][i] = u0_inter[i][j];
	r0[j][i] = r0[i][j];
	r02_inter[j][i] = r02_inter[i][j];
	rhoScl_inter[j][i] =rhoScl_inter[i][j];
	semp1_inter[j][i]=semp1_inter[i][j];
	offset[j][i] = offset[i][j];
	
	
	return cut[i][j];
}

/* ----------------------------------------------------------------------
 proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void PairSemp::write_restart(FILE *fp)
{
	std::cout<<"write_restart"<<std::endl;
	write_restart_settings(fp);
	
	int i,j;
//	for (i = 1; i <= atom->ntypes; i++)
//		for (j = i; j <= atom->ntypes; j++) {
//			fwrite(&setflag[i][j],sizeof(int),1,fp);
//			if (setflag[i][j]) {
//				fwrite(&d0[i][j],sizeof(double),1,fp);
//				fwrite(&alpha[i][j],sizeof(double),1,fp);
//				fwrite(&r0[i][j],sizeof(double),1,fp);
//				fwrite(&cut[i][j],sizeof(double),1,fp);
//			}
//		}
}

/* ----------------------------------------------------------------------
 proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void PairSemp::read_restart(FILE *fp)
{
	std::cout<<"read_restart"<<std::endl;
	read_restart_settings(fp);
	
	allocate();
	
	int i,j;
//	int me = comm->me;
//	for (i = 1; i <= atom->ntypes; i++)
//		for (j = i; j <= atom->ntypes; j++) {
//			if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
//			MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
//			if (setflag[i][j]) {
//				if (me == 0) {
//					fread(&d0[i][j],sizeof(double),1,fp);
//					fread(&alpha[i][j],sizeof(double),1,fp);
//					fread(&r0[i][j],sizeof(double),1,fp);
//					fread(&cut[i][j],sizeof(double),1,fp);
//				}
//				MPI_Bcast(&d0[i][j],1,MPI_DOUBLE,0,world);
//				MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
//				MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
//				MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
//			}
//		}
}

/* ----------------------------------------------------------------------
 proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void PairSemp::write_restart_settings(FILE *fp)
{
	std::cout<<"write_restart_settings"<<std::endl;	
//	fwrite(&cut_global,sizeof(double),1,fp);
//	fwrite(&offset_flag,sizeof(int),1,fp);
//	fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
 proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void PairSemp::read_restart_settings(FILE *fp)
{
	std::cout<<"read_restart_settings"<<std::endl;
//	if (comm->me == 0) {
//		fread(&cut_global,sizeof(double),1,fp);
//		fread(&offset_flag,sizeof(int),1,fp);
//		fread(&mix_flag,sizeof(int),1,fp);
//	}
//	MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
//	MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
//	MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairSemp::single(int i, int j, int itype, int jtype, int imol, int jmol, double rsq,
						 double factor_coul, double factor_lj,
						 double &fforce)
{
	std::cout<<"single"<<std::endl;
	double r,dr,dexp,phi;
	
	dr = rsq - r02[itype,jtype];
	if(imol == jmol){
		dexp = exp(-rhoScl[itype] * dr);
		fforce = factor_lj * semp1[itype]*(dexp*dexp - dexp);
		phi=u0[itype] * (dexp*dexp - 2.0*dexp) - offset[itype][jtype];
				
	}else{
		dexp = exp(-rhoScl_inter[itype][jtype]*dr);
		fforce = factor_lj * semp1_inter[itype][jtype]*(dexp*dexp - dexp);
		phi=u0_inter[itype][jtype] * (dexp*dexp - 2.0*dexp) - offset[itype][jtype];
		
	}
	
	return factor_lj*phi;
}
