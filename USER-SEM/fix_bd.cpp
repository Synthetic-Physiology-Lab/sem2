/*
 *  fix_bd.cpp
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 1/20/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 * 
 * ------------------------------
 * implements brownian dynamics
 *
 */

#include "stdio.h"
#include "string.h"

#include "stdlib.h"
#include "fix_bd.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBD::FixBD(LAMMPS *lmp, int narg, char **arg) :
Fix(lmp, narg, arg)
{
    if(narg != 4) error->all(FLERR," fix bd needs damping factor to be defined (same as for langevin");
    t_period = atof(arg[3]);
    dynamic_group_allow = 1;
    time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixBD::setmask()
{
    int mask = 0;
    mask |= INITIAL_INTEGRATE;
    mask |= FINAL_INTEGRATE;
    mask |= INITIAL_INTEGRATE_RESPA;
    mask |= FINAL_INTEGRATE_RESPA;
    mask |= PRE_FORCE;
    mask |= PRE_FORCE_RESPA;    
    return mask;
}

/* ---------------------------------------------------------------------- */


void FixBD::init(){
    dtv = update->dt;
    dtf = 0.5 * update->dt * force->ftm2v;
    
    if (strcmp(update->integrate_style,"respa") == 0)
	step_respa = ((Respa *) update->integrate)->step;
}
 
 /* ----------------------------------------------------------------------
    sets velocity to v0
    assumes it is called after langevin fix
 ------------------------------------------------------------------------- */

void FixBD::setup(int vflag){
    
    double dtfm;
    
    // update v of atoms in group
    double **v = atom->v;
    double **f = atom->f;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    if (igroup == atom->firstgroup) nlocal = atom->nfirst;
    
    if (rmass) {
	for (int i = 0; i < nlocal; i++)
	    if (mask[i] & groupbit) {
		dtfm = force->ftm2v * t_period / rmass[i];
		v[i][0] = dtfm * f[i][0];
		v[i][1] = dtfm * f[i][1];
		v[i][2] = dtfm * f[i][2];
	    }
    } else {
	for (int i = 0; i < nlocal; i++)
	    if (mask[i] & groupbit) {
		dtfm = force->ftm2v * t_period / mass[type[i]];
		v[i][0] = dtfm * f[i][0];
		v[i][1] = dtfm * f[i][1];
		v[i][2] = dtfm * f[i][2];
	    }
    }
}

/* ----------------------------------------------------------------------
 sets velocity to 0 in order to cancel dumping in langevin fix
 ------------------------------------------------------------------------- */

void FixBD::pre_force(int vflag)
{
    double **v = atom->v;
    double nlocal = atom->nlocal;
    for (int i=0; i < nlocal; i++){
	v[i][0]=0.0;
	v[i][1]=0.0;
	v[i][2]=0.0;
    }
}

/* ----------------------------------------------------------------------
 assume velocity updated at the end of timestep / before initial step
 bring x_i to x_i+1
 ------------------------------------------------------------------------- */

void FixBD::initial_integrate(int vflag)
{
    double dtfm;
    
    // update x of atoms in group
    
    double **x = atom->x;
    double **v = atom->v;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    if (igroup == atom->firstgroup) nlocal = atom->nfirst;
    
    for (int i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) {
	    x[i][0] += dtv * v[i][0];
	    x[i][1] += dtv * v[i][1];
	    x[i][2] += dtv * v[i][2];
	}  
}


/* ----------------------------------------------------------------------
    update velocity to v_i+1
 ------------------------------------------------------------------------- */

void FixBD::final_integrate()
{
    double dtfm;
    
    // update v of atoms in group
    double **v = atom->v;
    double **f = atom->f;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    if (igroup == atom->firstgroup) nlocal = atom->nfirst;
    
    if (rmass) {
	for (int i = 0; i < nlocal; i++)
	    if (mask[i] & groupbit) {
		dtfm = force->ftm2v * t_period / rmass[i];
		v[i][0] = dtfm * f[i][0];
		v[i][1] = dtfm * f[i][1];
		v[i][2] = dtfm * f[i][2];
	    }
	
    } else {
	for (int i = 0; i < nlocal; i++)
	    if (mask[i] & groupbit) {
		dtfm = force->ftm2v * t_period / mass[type[i]];
		v[i][0] = dtfm * f[i][0];
		v[i][1] = dtfm * f[i][1];
		v[i][2] = dtfm * f[i][2];
	    }
    }
    
}

/* ---------------------------------------------------------------------- */

void FixBD::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
    dtv = step_respa[ilevel];
    dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
    
    // innermost level - NVE update of v and x
    // all other levels - NVE update of v
    
    if (ilevel == 0) initial_integrate(vflag);
    else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixBD::final_integrate_respa(int ilevel, int iloop)
{
    dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
    final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixBD::reset_dt()
{
    dtv = update->dt;
    dtf = 0.5 * update->dt * force->ftm2v;
}
