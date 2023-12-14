/*
 *  dump_xyzCell.h
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 1/12/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 *  ++++++++++
 *  based on dump_xyz
 *  ++++++++++
 */

#ifdef DUMP_CLASS

DumpStyle(xyzCell,DumpXYZCell)

#else

#ifndef LMP_DUMP_XYZCELL_H
#define LMP_DUMP_XYZCELL_H

#include "dump.h"

namespace LAMMPS_NS {
    
    class DumpXYZCell : public Dump {
    public:
	DumpXYZCell(class LAMMPS *, int, char**);
	virtual ~DumpXYZCell() {}
	
    protected:
	void init_style();
	void write_header(bigint);
	int count();
	void pack(tagint *);
	void write_data(int, double *);
    };
    
}

#endif
#endif