/*
 *  dump_xyzSemEM.h
 *  LAMMPS_xcode
 *
 *  Created by Florian Milde on 1/12/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 *  ++++++++++
 *  based on dump_xyzCell
 *  ++++++++++
 */

#ifdef DUMP_CLASS

DumpStyle(xyzSemEM,DumpXYZSemEM)

#else

#ifndef LMP_DUMP_XYZSEMEM_H
#define LMP_DUMP_XYZSEMEM_H

#include "dump.h"

namespace LAMMPS_NS {
    
    class DumpXYZSemEM : public Dump {
    public:
	DumpXYZSemEM(class LAMMPS *, int, char**);
	virtual ~DumpXYZSemEM() {}
	
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