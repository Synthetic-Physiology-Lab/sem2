/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
 
   ++++++ 
   + Adapted atom style for SEM simulations
   + Based on atom style full
   + Florian Milde, mildef@ethz.ch
   ++++++  
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS

AtomStyle(sem,AtomVecSem)

#else

#ifndef LMP_ATOM_VEC_SEM_H
#define LMP_ATOM_VEC_SEM_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecSem : public AtomVec {
 public:
  AtomVecSem(class LAMMPS *);
  virtual ~AtomVecSem() {}
  virtual void grow(int);
  virtual void grow_reset();
  virtual void copy(int, int, int);
  virtual int pack_comm(int, int *, double *, int, int *);
  virtual int pack_comm_vel(int, int *, double *, int, int *);
  virtual void unpack_comm(int, int, double *);
  virtual void unpack_comm_vel(int, int, double *);
  virtual int pack_reverse(int, int, double *);
  virtual void unpack_reverse(int, int *, double *);
  virtual int pack_border(int, int *, double *, int, int *);
  virtual int pack_border_vel(int, int *, double *, int, int *);
  virtual int pack_border_hybrid(int, int *, double *);  //@@@???
  virtual void unpack_border(int, int, double *);
  virtual void unpack_border_vel(int, int, double *);
  virtual int unpack_border_hybrid(int, int, double *);  //@@@???
  virtual int pack_exchange(int, double *);
  virtual int unpack_exchange(double *);
  virtual int size_restart();
  virtual int pack_restart(int, double *);
  virtual int unpack_restart(double *);
  virtual void create_atom(int, double *);
  virtual void data_atom(double *, imageint, char **);
  virtual int data_atom_hybrid(int, char **);//@@@???
  virtual void pack_data(double **);
  virtual int pack_data_hybrid(int, double *);
  virtual void write_data(FILE *, int, double **);
  virtual int write_data_hybrid(FILE *, double *);
  virtual double memory_usage();

 protected:
  tagint *tag;
  int *type,*mask;
  imageint *image;
  double **x,**v,**f;
  double *q;
  int **nspecial;
  tagint **special;//@@@???
  int *num_bond;
  int **bond_type;
  tagint **bond_atom;
  //@@@SEM
  double *p;   // p holds polymerization factor
  int *cell;      // cell id (much like molecule id
  int *cellpart;  // 0=cytoplasm; 1=membrane; 2=nucleus	
};

}

#endif
#endif
