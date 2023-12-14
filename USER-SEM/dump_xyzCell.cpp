/*
 *  dump_xyzCell.cpp
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 1/12/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 *  ++++++++++
 *  based on dump_xyz
 *  ++++++++++ 
 */

#include "string.h"
#include "dump_xyzCell.h"
#include "atom.h"
#include "group.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpXYZCell::DumpXYZCell(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
    if (narg != 5) error->all(FLERR,"Illegal dump xyz command");
    if (binary || multiproc) error->all(FLERR,"Invalid dump xyz filename");
    if (!atom->sem_flag) error->all(FLERR,"dump style only supported for sem package");
    
    size_one = 5;
    sort_flag = 1;
    sortcol = 0;
    
    char *str = (char *) "%d %g %g %g";
    int n = strlen(str) + 1;
    format_default = new char[n];
    strcpy(format_default,str);
}

/* ---------------------------------------------------------------------- */

void DumpXYZCell::init_style()
{
    delete [] format;
    char *str;
    if (format_line_user) str = format_line_user;
    else str = format_default;
    
    int n = strlen(str) + 2;
    format = new char[n];
    strcpy(format,str);
    strcat(format,"\n");
    
    // open single file, one time only
    
    if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

void DumpXYZCell::write_header(bigint n)
{
    if (me == 0) {
	fprintf(fp,BIGINT_FORMAT "\n",n);
	fprintf(fp,"Atoms\n");
    }
}

/* ---------------------------------------------------------------------- */

int DumpXYZCell::count()
{
    if (igroup == 0) return atom->nlocal;
    
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    
    int m = 0;
    for (int i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) m++;
    return m;
}

/* ---------------------------------------------------------------------- */

void DumpXYZCell::pack(tagint *ids)
{
    int m,n;
    
    tagint *tag = atom->tag;
    int *type = atom->type;
    int *cell = atom->cell;
    int *mask = atom->mask;
    double **x = atom->x;
    int nlocal = atom->nlocal;
    
    m = n = 0;
    for (int i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) {
	    buf[m++] = tag[i];
	    buf[m++] = cell[i];
	    buf[m++] = x[i][0];
	    buf[m++] = x[i][1];
	    buf[m++] = x[i][2];
	    if (ids) ids[n++] = tag[i];
	}
}

/* ---------------------------------------------------------------------- */

void DumpXYZCell::write_data(int n, double *mybuf)
{
    int m = 0;
    for (int i = 0; i < n; i++) {
	fprintf(fp,format,
		static_cast<int> (mybuf[m+1]),mybuf[m+2],mybuf[m+3],mybuf[m+4]);
	m += size_one;
    }
}
