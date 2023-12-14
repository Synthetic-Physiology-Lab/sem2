/*
 *  pari_semp.cpp
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 12/9/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_semextp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

#include "update.h"
using namespace LAMMPS_NS;

#define EPS 1e-5

// define __SEM_DEBUG__ to get additional printfs
//#define __SEM_DEBUG__

/* ---------------------------------------------------------------------- */

PairSemextp::PairSemextp(LAMMPS *lmp) : Pair(lmp) {}

/* ----------------------------------------------------------------------
// what was allocated

 // what to allocate
 
 // genral
 setflag[ii,jj]
 cutsq[ii,jj]
 
 r0[cellType];
 
 // intra

 u0[cellType]
 rho[cellType]
 alpha[cellType]
 rhoScl[cellType]
 semp1[cellType]
 cut_loc[cellType] 
 
 // inter tmp
 u0_inter[ii,jj]
 rho_inter[ii,jj]
 alpha_inter[ii,jj]
 rhoScl_inter[ii,jj]
 semp1_inter[ii,jj]
 cut[ii,jj]

--------------------------------------------------------------------------*/

PairSemextp::~PairSemextp()
{
    if (allocated) {
        // general
        memory->destroy(setflag);
        memory->destroy(cutsq);
        
        memory->destroy(r0);
        
        // intra
        memory->destroy(u0);
        memory->destroy(rho);
        memory->destroy(alpha);
        memory->destroy(rhoScl);
        memory->destroy(semp1);
        memory->destroy(cut_loc);
        
        // inter
        memory->destroy(u0_inter);
        memory->destroy(rho_inter);
        memory->destroy(alpha_inter);       
        memory->destroy(rhoScl_inter);
        memory->destroy(semp1_inter);
        memory->destroy(cut);
        
        // ??
        memory->destroy(offset);
    }
}

/* ---------------------------------------------------------------------- */

void PairSemextp::compute(int eflag, int vflag)
{
    int i,j,ii,jj,inum,jnum,itype,jtype,icell,jcell;
    double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair,ip,jp;
    double rsq,pFac,E,factor_lj;
    int *ilist,*jlist,*numneigh,**firstneigh;
    
    evdwl = 0.0;
    if (eflag || vflag) ev_setup(eflag,vflag);
    else evflag = vflag_fdotr = 0;
    
    double **x = atom->x;
    double **f = atom->f;
    int *type = atom->type;
    int *cell  = atom->cell;
    double *p     = atom->p;
    int nlocal = atom->nlocal;
    int nall = nlocal + atom->nghost;
    double *special_lj = force->special_lj;
    int newton_pair = force->newton_pair;
    
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    
    // loop over neighbors of my atoms
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
#ifdef __SEM_DEBUG__
        printf("++i" BIGINT_FORMAT " | proc: %d, particle %d, location: %f %f %f\n",update->ntimestep, me, atom->tag[i],x[i][0],x[i][1],x[i][2]);
#endif


        jlist = firstneigh[i];
        jnum = numneigh[i];
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
            rsq = (delx*delx + dely*dely + delz*delz);
            if (rsq < cutsq[itype][jtype]) {

                pFac = MIN(ip,jp);
                rsq/=(pFac*pFac);

                //double p=MIN(ip,jp);
                //if(p<=0.0) continue;
                //rsq/=(p*p);

                // inter
                if(icell == jcell){
                    E = exp(rho[itype]-rhoScl[itype]*rsq);

                    fpair = factor_lj * semp1[itype]*(2.*E*E-alpha[itype]*E)*pFac;

                    f[i][0] += delx*fpair;
                    f[i][1] += dely*fpair;
                    f[i][2] += delz*fpair;
                    if (newton_pair || j < nlocal) {
                        f[j][0] -= delx*fpair;
                        f[j][1] -= dely*fpair;
                        f[j][2] -= delz*fpair;
                    }
                    
                    if (eflag) {
                        evdwl = u0[itype] * (E*E - alpha[itype]*E) -
                        offset[itype][jtype];
                        evdwl *= factor_lj;
                    }
                    
                    
                }else{
                    E = exp(rho_inter[itype][jtype]-rhoScl_inter[itype][jtype]*rsq);
                    fpair = factor_lj * semp1_inter[itype][jtype]*(2.*E*E-alpha_inter[itype][jtype]*E)*pFac;
                    
                    f[i][0] += delx*fpair;
                    f[i][1] += dely*fpair;
                    f[i][2] += delz*fpair;
                    if (newton_pair || j < nlocal) {
                        f[j][0] -= delx*fpair;
                        f[j][1] -= dely*fpair;
                        f[j][2] -= delz*fpair;
                    }
                    if (eflag) {
                        evdwl = u0_inter[itype][jtype] * (E*E - alpha_inter[itype][jtype]) -
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
}

/* ----------------------------------------------------------------------
 allocate all arrays 
 
 // what to allocate
 
 // genral
 setflag[ii,jj]
 cutsq[ii,jj]
 
 r0[cellType];
 
 // intra
 
 u0[cellType]
 rho[cellType]
 alpha[cellType]
 rhoScl[cellType]
 semp1[cellType]
 cut_loc[cellType] 
 
 // inter tmp
 u0_inter[ii,jj]
 rho_inter[ii,jj]
 alpha_inter[ii,jj]
 rhoScl_inter[ii,jj]
 semp1_inter[ii,jj]
 cut[ii,jj]
 ------------------------------------------------------------------------- */

void PairSemextp::allocate()
{
    allocated = 1;
    int n = atom->ntypes; // cell types
    
    // general
    memory->create(setflag,n+1,n+1,"pair:setflag");
    for (int i = 1; i <= n; i++)
        for (int j = i; j <= n; j++)
            setflag[i][j] = 0;
    
    memory->create(cutsq,n+1,n+1,"pair:cutsq");
    memory->create(r0,n,"pair:r0");
    
    //intra
    memory->create(rho,n,"pair:rho");
    memory->create(u0,n,"pair:u0");
    memory->create(alpha,n,"pair:alpha");
    memory->create(rhoScl,n,"pair:rhoScl");
    memory->create(semp1,n,"pair:semp1");
    memory->create(cut_loc,n,"pair:cut_loc");
    
    //intra tmp
    memory->create(rho_inter,n+1,n+1,"pair:rho_inter");
    memory->create(u0_inter,n+1,n+1,"pair:u0_inter");
    memory->create(alpha_inter,n+1,n+1,"pair:alpha_inter");
    memory->create(rhoScl_inter,n+1,n+1,"pair:rhoScl_inter");
    memory->create(semp1_inter,n+1,n+1,"pair:semp1_inter");
    memory->create(cut,n+1,n+1,"pair:cut");
    
    //??
    memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
 global settings 
 ------------------------------------------------------------------------- */

void PairSemextp::settings(int narg, char **arg)
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
    -0) atom_type   (CELL TYPE)
    -1) place holder / not considered
    -2) u0_intra    (u_0)
    -3) u0_inter
    -4) rho_intra   [internally, store alpha=rho/req2]
    -5) rho_inter
    -6) alpha_intra
    -7) alpha_inter
    -8) r0          (r_eq)[store r0^2]
    -9) cutoff      (optional)
we store internally
    rho(_inter):    rho
    u0(_inter):     u0
    alpha(_inter):  alpha
    rhoScl(_inter): rho/r0^2
    semp1(_inter):  2*u0*rho/r0^2
    cut_loc/cut:    cutoff
 
    r0:             r_eq

 ------------------------------------------------------------------------- */

void PairSemextp::coeff(int narg, char **arg)
{
    int ii,jj;
    if (narg < 9 || narg > 10) error->all(FLERR,"Incorrect args for pair coefficients");
    if (!allocated) allocate();
    
    int cellType=utils::inumeric(FLERR,arg[0],false,lmp);
    int spare = utils::inumeric(FLERR,arg[1],false,lmp);
    printf("PAIRCOEFF_SEMEXT: second parameter %d ignored.. mixing rules applied\n",spare);
    
    double u0_intra_one = utils::numeric(FLERR,arg[2],false,lmp);
    double u0_inter_one = utils::numeric(FLERR,arg[3],false,lmp);
    double rho_intra_one = utils::numeric(FLERR,arg[4],false,lmp);
    double rho_inter_one = utils::numeric(FLERR,arg[5],false,lmp);
    double alpha_intra_one = utils::numeric(FLERR,arg[6],false,lmp);
    double alpha_inter_one = utils::numeric(FLERR,arg[7],false,lmp);    
    double r0_one = utils::numeric(FLERR,arg[8],false,lmp);
       
    double cut_one = cut_global;
    if (narg == 10) cut_one = utils::numeric(FLERR,arg[9],false,lmp);

    // set general coefficients
    setflag[cellType][cellType]=1;
    r0[cellType]=r0_one;

    
    // set intra coefficients   
    u0[cellType]=u0_intra_one;
    rho[cellType]=rho_intra_one;
    alpha[cellType]=alpha_intra_one;
    cut_loc[cellType]=cut_one;
    
    // store inter coefficients
    u0_inter[cellType][cellType]=u0_inter_one;
    rho_inter[cellType][cellType]=rho_inter_one;
    alpha_inter[cellType][cellType]=alpha_inter_one;
    cut[cellType][cellType]=cut_one;
    
    printf(" set coeffs: u0 %f, u0_inter %f, rho0 %f, rho0_inter %f, alpha %f, alpha_inter %f, r0 %f\n",
           u0[cellType],u0_inter[cellType][cellType],rho[cellType],rho_inter[cellType][cellType],alpha[cellType],alpha_inter[cellType][cellType],r0[cellType]);

/*
    // calculated coefficients
    double r02=r0_one*r0_one;
    // intra
    rhoScl[cellType]=rho_intra_one/r02;
    semp1[cellType]=2.*rho_intra_one*u0_intra_one/r02;

    // inter
    rhoScl_inter[cellType][cellType]=rho_inter_one/r02;
    semp1_inter[cellType][cellType]=2.0*u0_inter_one*rhoScl_inter[cellType][cellType];
*/

    // update intra coefficient table
    for (int i=1;i<=atom->ntypes;i++) {
        // only if other potential is set
        if(setflag[i][i]) {
            setflag[cellType][i]=1;
/*
            ii=MIN(i,cellType);
            jj=MAX(i,cellType);
            u0_inter[ii][jj]=0.5*(u0_inter[ii][ii]+u0_inter[jj][jj]);
            rho_inter[ii][jj]=0.5*(rho_inter[ii][ii]+u0_inter[jj][jj]);
            alpha_inter[ii][jj]=0.5*(alpha_inter[ii][ii]+alpha_inter[jj][jj]);          
            rhoScl_inter[ii][jj]=4.0*rho_inter[ii][jj]/((r0[ii]+r0[jj])*(r0[ii]+r0[jj]));
            semp1_inter[ii][jj]=2.0*u0_inter[ii][jj]*rhoScl_inter[ii][jj];
            cut[ii][jj]=MAX(cut_loc[ii],cut_loc[jj]);
*/
        }
    };
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSemextp::init_one(int i, int j)
{      
    printf("init_one for %d and %d\n",i,j);
    double r02;
    if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
    if(i==j){
        // also set intra derived parameters
        r02=r0[i]*r0[i];
        rhoScl[i]=rho[i]/r02;
        semp1[i]=2.*rho[i]*u0[i]/r02;
    }
        
            
    offset[i][j] = 0.0;
    // scaled
    if(i!=j){
        u0_inter[i][j]=0.5*(u0_inter[i][i]+u0_inter[j][j]);
        rho_inter[i][j]=0.5*(rho_inter[i][i]+rho_inter[j][j]);
        alpha_inter[i][j]=0.5*(alpha_inter[i][i]+alpha_inter[j][j]);
        cut[i][j]=MAX(cut_loc[i],cut_loc[j]);
    }
    // derived
    rhoScl_inter[i][j]=4.0*rho_inter[i][j]/((r0[i]+r0[j])*(r0[i]+r0[j]));
    semp1_inter[i][j]=2.0*u0_inter[i][j]*rhoScl_inter[i][j];
    
    // copy to j i
    if(i!=j){
        offset[j][i]=offset[i][j];
        rho_inter[j][i]=rho_inter[i][j];
        alpha_inter[j][i]=alpha_inter[i][j];
        rhoScl_inter[j][i]=rhoScl_inter[i][j];
        semp1_inter[j][i]=semp1_inter[i][j];
        cut[j][i]=cut[i][j];
    }
    return cut[i][j];
}

/* ----------------------------------------------------------------------
 proc 0 writes to restart file
 ------------------------------------------------------------------------- */
//@@@SEM RESTART NOT TESTED
void PairSemextp::write_restart(FILE *fp)
{
    write_restart_settings(fp);
    
    int i,j;
    for (i = 1; i <= atom->ntypes; i++){
        fwrite(&setflag[i][i],sizeof(int),1,fp);
        if(setflag[i][i]){
            //intra
            fwrite(&u0[i],sizeof(double),1,fp);
            fwrite(&rho[i],sizeof(double),1,fp);
            fwrite(&alpha[i],sizeof(double),1,fp);
            fwrite(&cut_loc[i],sizeof(double),1,fp);
            //inter
            fwrite(&u0_inter[i][i],sizeof(double),1,fp);
            fwrite(&rho_inter[i][i],sizeof(double),1,fp);
            fwrite(&alpha_inter[i][i],sizeof(double),1,fp);
            fwrite(&cut[i][i],sizeof(double),1,fp);
            
        }
        for (j = i+1; j <= atom->ntypes; j++) {
            fwrite(&setflag[i][j],sizeof(int),1,fp);
        }
    }
}

/* ----------------------------------------------------------------------
 proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */
//@@@SEM RESTART NOT TESTED
void PairSemextp::read_restart(FILE *fp)
{
    read_restart_settings(fp);
    
    allocate();
    
    int i,j;
    int me = comm->me;
    for (i = 1; i <= atom->ntypes; i++){
        if (me == 0) fread(&setflag[i][i],sizeof(int),1,fp);
        MPI_Bcast(&setflag[i][i],1,MPI_INT,0,world);
        if (setflag[i][i]) {
            if (me == 0){
                //intra
                fread(&u0[i],sizeof(double),1,fp);
                fread(&rho[i],sizeof(double),1,fp);
                fread(&alpha[i],sizeof(double),1,fp);
                fread(&cut_loc[i],sizeof(double),1,fp);
                //inter
                fread(&u0_inter[i][i],sizeof(double),1,fp);
                fread(&rho_inter[i][i],sizeof(double),1,fp);
                fread(&alpha_inter[i][i],sizeof(double),1,fp);
                fread(&cut[i][i],sizeof(double),1,fp);
            }
            //intra
            MPI_Bcast(&u0[i],1,MPI_DOUBLE,0,world);
            MPI_Bcast(&rho[i],1,MPI_DOUBLE,0,world);
            MPI_Bcast(&alpha[i],1,MPI_DOUBLE,0,world);
            MPI_Bcast(&cut_loc[i],1,MPI_DOUBLE,0,world);
            //inter
            MPI_Bcast(&u0_inter[i][i],1,MPI_DOUBLE,0,world);
            MPI_Bcast(&rho_inter[i][i],1,MPI_DOUBLE,0,world);
            MPI_Bcast(&alpha_inter[i][i],1,MPI_DOUBLE,0,world);
            MPI_Bcast(&cut[i][i],1,MPI_DOUBLE,0,world);
        }
        for (j = i+1; j <= atom->ntypes; j++) {
            if(me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
            MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
        }
                
    }
}

/* ----------------------------------------------------------------------
 proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void PairSemextp::write_restart_settings(FILE *fp)
{   
    fwrite(&cut_global,sizeof(double),1,fp);
    fwrite(&offset_flag,sizeof(int),1,fp);
    fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
 proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void PairSemextp::read_restart_settings(FILE *fp)
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
double PairSemextp::single(int i, int j, int itype, int jtype, int icell, int jcell,double poly, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
    double r,E,phi;
            
    ////////
    r=sqrt(rsq);
    rsq/=(poly*poly+EPS);
            
    if(icell == jcell){
        E = exp(rho[itype]-rhoScl[itype]*rsq);
        fforce = factor_lj * semp1[itype]*(2.*E*E - alpha[itype]*E)/r;
        phi=u0[itype] * (E*E - alpha[itype]*E) - offset[itype][jtype];
        
    }else{
        E = exp(rho_inter[itype][jtype]-rhoScl_inter[itype][jtype]*rsq);
        fforce = factor_lj * semp1_inter[itype][jtype]*(2.*E*E - alpha_inter[itype][jtype]*E)/r;
        phi=u0_inter[itype][jtype] * (E*E - alpha_inter[itype][jtype]*E) - offset[itype][jtype];
    }
    return factor_lj*phi;
}
