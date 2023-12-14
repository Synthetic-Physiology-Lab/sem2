/*
 *  atom_vec_semEM.cpp
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 1/26/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 *  ++++++++++++++++++++++++++++++++++++++
 *   - to use with pair_semEM
 *   - based on atom_vec_sem
 *  ++++++++++++++++++++++++++++++++++++++
 */

#include "stdlib.h"
#include "atom_vec_semEM.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecSemEM::AtomVecSemEM(LAMMPS *lmp) : AtomVecSem(lmp)
{
  size_border = 11;   //ok
  size_data_atom = 10; //ok
  xcol_data = 8;      //ok
  atom->semDens_flag = 1;
}

/* ----------------------------------------------------------------------
 grow atom arrays
 n = 0 grows arrays by a chunk
 n > 0 allocates arrays to size n 
 ------------------------------------------------------------------------- */

void AtomVecSemEM::grow(int n)
{
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");

  q = memory->grow(atom->q,nmax,"atom:q");
  //molecule = memory->grow(atom->molecule,nmax,"atom:molecule");

  nspecial = memory->grow(atom->nspecial,nmax,3,"atom:nspecial");
  special = memory->grow(atom->special,nmax,atom->maxspecial,"atom:special");

  num_bond = memory->grow(atom->num_bond,nmax,"atom:num_bond");
  bond_type = memory->grow(atom->bond_type,nmax,atom->bond_per_atom,
			   "atom:bond_type");
  bond_atom = memory->grow(atom->bond_atom,nmax,atom->bond_per_atom,
			   "atom:bond_atom");

  //@@@SEM
  p = memory->grow(atom->p,nmax,"atom:p");
  cell = memory->grow(atom->cell,nmax,"atom:cell");
  cellpart = memory->grow(atom->cellpart,nmax,"atom:cellpart");
  dens = atom->dens = memory->grow(atom->dens,nmax,"atom:dens");
  
  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) 
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
 reset local array ptrs
 ------------------------------------------------------------------------- */

void AtomVecSemEM::grow_reset()
{
  AtomVecSem::grow_reset();
  dens = atom->dens;
}

/* ---------------------------------------------------------------------- */

void AtomVecSemEM::copy(int i, int j, int delflag)
{
  int k;

  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];
	
  q[j] = q[i];

  num_bond[j] = num_bond[i];
  for (k = 0; k < num_bond[j]; k++) {
    bond_type[j][k] = bond_type[i][k];
    bond_atom[j][k] = bond_atom[i][k];
  }

  nspecial[j][0] = nspecial[i][0];
  nspecial[j][1] = nspecial[i][1];
  nspecial[j][2] = nspecial[i][2];
  for (k = 0; k < nspecial[j][2]; k++) special[j][k] = special[i][k];

  //@@@SEM
  p[j] = p[i];
  cell[j] = cell[i];
  cellpart[j] = cellpart[i];
  dens[j] = dens[i];
	
  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) 
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecSemEM::pack_border(int n, int *list, double *buf,
			    int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = q[j];
      buf[m++] = p[j];
      buf[m++] = ubuf(cell[j]).d;
      buf[m++] = ubuf(cellpart[j]).d;
      buf[m++] = dens[j];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = q[j];
      buf[m++] = p[j];
      buf[m++] = ubuf(cell[j]).d;
      buf[m++] = ubuf(cellpart[j]).d;
      buf[m++] = dens[j];
    }
  }
  
  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);
      
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSemEM::pack_border_vel(int n, int *list, double *buf,
				 int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = q[j];
      buf[m++] = p[j];
      buf[m++] = ubuf(cell[j]).d;
      buf[m++] = ubuf(cellpart[j]).d;
      buf[m++] = dens[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = q[j];
        buf[m++] = p[j];
        buf[m++] = ubuf(cell[j]).d;
        buf[m++] = ubuf(cellpart[j]).d;
        buf[m++] = dens[j];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = q[j];
        buf[m++] = p[j];
        buf[m++] = ubuf(cell[j]).d;
        buf[m++] = ubuf(cellpart[j]).d;
        buf[m++] = dens[j];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSemEM::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = q[j];
    buf[m++] = p[j];
    buf[m++] = ubuf(cell[j]).d;
    buf[m++] = ubuf(cellpart[j]).d;
    buf[m++] = dens[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSemEM::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    q[i] = buf[m++];
    p[i] = buf[m++];
    cell[i] = (int) ubuf(buf[m++]).i;
    cellpart[i] = (int) ubuf(buf[m++]).i;
    dens[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecSemEM::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    q[i] = buf[m++];
    p[i] = buf[m++];
    cell[i] = (int) ubuf(buf[m++]).i;
    cellpart[i] = (int) ubuf(buf[m++]).i;
    dens[i] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecSemEM::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    q[i] = buf[m++];
    p[i] = buf[m++];
    cell[i] = (int) ubuf(buf[m++]).i;
    cellpart[i] = (int) ubuf(buf[m++]).i;
    dens[i] = buf[m++];
  }
  return m;
}

/* ----------------------------------------------------------------------
 pack data for atom I for sending to another proc
 xyz must be 1st 3 values, so comm::exchange() can test on them 
 ------------------------------------------------------------------------- */

int AtomVecSemEM::pack_exchange(int i, double *buf)
{
  int k;

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;

  buf[m++] = q[i];
  buf[m++] = p[i];
  buf[m++] = ubuf(cell[i]).d;
  buf[m++] = ubuf(cellpart[i]).d;
  buf[m++] = dens[i];

  buf[m++] = ubuf(num_bond[i]).d;
  for (k = 0; k < num_bond[i]; k++) {
    buf[m++] = ubuf(bond_type[i][k]).d;
    buf[m++] = ubuf(bond_atom[i][k]).d;
  }

  buf[m++] = ubuf(nspecial[i][0]).d;
  buf[m++] = ubuf(nspecial[i][1]).d;
  buf[m++] = ubuf(nspecial[i][2]).d;
  for (k = 0; k < nspecial[i][2]; k++) buf[m++] = ubuf(special[i][k]).d;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSemEM::unpack_exchange(double *buf)
{
  int k;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;

  q[nlocal] = buf[m++];
  p[nlocal] = buf[m++];
  cell[nlocal] = (int) ubuf(buf[m++]).i;
  cellpart[nlocal] = (int) ubuf(buf[m++]).i;
  dens[nlocal] = buf[m++];

  num_bond[nlocal] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < num_bond[nlocal]; k++) {
    bond_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    bond_atom[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }

  nspecial[nlocal][0] = (int) ubuf(buf[m++]).i;
  nspecial[nlocal][1] = (int) ubuf(buf[m++]).i;
  nspecial[nlocal][2] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < nspecial[nlocal][2]; k++)
    special[nlocal][k] = (tagint) ubuf(buf[m++]).i;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) 
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
 size of restart data for all atoms owned by this proc
 include extra data stored by fixes
 ------------------------------------------------------------------------- */

int AtomVecSemEM::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 0;
  for (i = 0; i < nlocal; i++)
    n += 17 + 2*num_bond[i];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++) 
      for (i = 0; i < nlocal; i++)
	    n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
 pack atom I's data for restart file including extra quantities
 xyz must be 1st 3 values, so that read_restart can test on them
 molecular types may be negative, but write as positive   
 ------------------------------------------------------------------------- */

int AtomVecSemEM::pack_restart(int i, double *buf)
{
  int k;

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  buf[m++] = q[i];
  buf[m++] = p[i];
  buf[m++] = ubuf(cell[i]).d;
  buf[m++] = ubuf(cellpart[i]).d;
  buf[m++] = dens[i];
  
  buf[m++] = ubuf(num_bond[i]).d;
  for (k = 0; k < num_bond[i]; k++) {
    buf[m++] = ubuf(MAX(bond_type[i][k],-bond_type[i][k])).d;
    buf[m++] = ubuf(bond_atom[i][k]).d;
  }
  
  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
 unpack data for one atom from restart file including extra quantities
 ------------------------------------------------------------------------- */

int AtomVecSemEM::unpack_restart(double *buf)
{
  int k;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  q[nlocal] = buf[m++];
  p[nlocal] = buf[m++];
  cell[nlocal] = (int) ubuf(buf[m++]).i;
  cellpart[nlocal] = (int) ubuf(buf[m++]).i;
  dens[nlocal] = buf[m++];
    
  num_bond[nlocal] = (int) ubuf(buf[m++]).i;
  for (k = 0; k < num_bond[nlocal]; k++) {
    bond_type[nlocal][k] = (int) ubuf(buf[m++]).i;
    bond_atom[nlocal][k] = (tagint) ubuf(buf[m++]).i;
  }
  
  nspecial[nlocal][0] = nspecial[nlocal][1] = nspecial[nlocal][2] = 0;

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
 create one atom of itype at coord with p=1
 set other values to defaults
 ------------------------------------------------------------------------- */

void AtomVecSemEM::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
    ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  q[nlocal] = 0.0;
  p[nlocal] = 1.0;
  cell[nlocal] = 0;
  cellpart[nlocal] = 0;
  dens[nlocal] = -1.0;
	
  num_bond[nlocal] = 0;
  nspecial[nlocal][0] = nspecial[nlocal][1] = nspecial[nlocal][2] = 0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
 
   FORMAT: atomID, cellID, atomType, q, p, cellpart, dens, x, y, z
   note: http://lammps.sandia.gov/doc/read_data.html
------------------------------------------------------------------------- */
//@@@???
void AtomVecSemEM::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);

  cell[nlocal] = atoi(values[1]);

  type[nlocal] = atoi(values[2]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  q[nlocal] = atof(values[3]);
  p[nlocal] = atof(values[4]);
  cellpart[nlocal] = atoi(values[5]);
  dens[nlocal] = atof(values[6]);

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  num_bond[nlocal] = 0;
	
  atom->nlocal++;
}

/* ----------------------------------------------------------------------
 unpack hybrid quantities from one line in Atoms section of data file
 initialize other atom quantities for this sub-style
 ------------------------------------------------------------------------- */
//@@@???
int AtomVecSemEM::data_atom_hybrid(int nlocal, char **values)
{
  cell[nlocal] = atoi(values[1]);
  cellpart[nlocal] = atoi(values[2]);
  q[nlocal] = atof(values[3]);
  dens[nlocal] = atof(values[4]);

  num_bond[nlocal] = 0;
  return 2;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecSemEM::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(cell[i]).d;
    buf[i][2] = ubuf(type[i]).d;
    buf[i][3] = q[i];
    buf[i][4] = x[i][0];
    buf[i][5] = x[i][1];
    buf[i][6] = x[i][2];
    buf[i][7] = ubuf(cellpart[i]).d;
    buf[i][8] = p[i];
    buf[i][9] = dens[i];
    buf[i][10] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][11] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][12] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecSemEM::pack_data_hybrid(int i, double *buf)
{
  buf[0] = ubuf(cell[i]).d;
  buf[1] = q[i];
  buf[2] = ubuf(cellpart[i]).d;
  buf[3] = p[i];
  buf[4] = dens[i];
  return 5;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecSemEM::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %d %d %-1.16e %-1.16e %-1.16e %-1.16e %d %-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            (int) ubuf(buf[i][2]).i,
            buf[i][3],buf[i][4],buf[i][5],buf[i][6],
            (int) ubuf(buf[i][7]).i,buf[i][8],buf[i][9],
            (int) ubuf(buf[i][10]).i,(int) ubuf(buf[i][11]).i,
            (int) ubuf(buf[i][12]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecSemEM::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %d %-1.16e %d %-1.16e %-1.16e",
          (int) ubuf(buf[0]).i, buf[1], (int) ubuf(buf[2]).i, buf[3], buf[4]);
  return 5;
}

/* ----------------------------------------------------------------------
 return # of bytes of allocated memory 
 ------------------------------------------------------------------------- */
double AtomVecSemEM::memory_usage()
{
  double bytes = 0.0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("q")) bytes += memory->usage(q,nmax);
  
  if (atom->memcheck("p")) bytes += memory->usage(p,nmax);
  if (atom->memcheck("cell")) bytes += memory->usage(cell,nmax);
  if (atom->memcheck("cellpart")) bytes += memory->usage(cellpart,nmax);
  if (atom->memcheck("dens")) bytes += memory->usage(dens,nmax);
  
  if (atom->memcheck("nspecial")) bytes += memory->usage(nspecial,nmax,3);
  if (atom->memcheck("special")) 
    bytes += memory->usage(special,nmax,atom->maxspecial);

  if (atom->memcheck("num_bond")) bytes += memory->usage(num_bond,nmax);
  if (atom->memcheck("bond_type")) 
    bytes += memory->usage(bond_type,nmax,atom->bond_per_atom);
  if (atom->memcheck("bond_atom")) 
    bytes += memory->usage(bond_atom,nmax,atom->bond_per_atom);

  return bytes;
}
